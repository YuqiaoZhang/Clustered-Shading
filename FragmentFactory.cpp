// Copyright 2014 Intel Corporation
// All Rights Reserved
//
// Permission is granted to use, copy, distribute and prepare derivative works of this
// software for any purpose and without fee, provided, that the above copyright notice
// and this statement appear in all copies.  Intel makes no representations about the
// suitability of this software for any purpose.  THIS SOFTWARE IS PROVIDED "AS IS."
// INTEL SPECIFICALLY DISCLAIMS ALL WARRANTIES, EXPRESS OR IMPLIED, AND ALL LIABILITY,
// INCLUDING CONSEQUENTIAL AND OTHER INDIRECT DAMAGES, FOR THE USE OF THIS SOFTWARE,
// INCLUDING LIABILITY FOR INFRINGEMENT OF ANY PROPRIETARY RIGHTS, AND INCLUDING THE
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  Intel does not
// assume any responsibility for any errors which may appear in this software nor any
// responsibility to update it.

#include "FragmentFactory.h"
#include <algorithm>

#include "ISPC/kernels.h"

class FragmentFactory
{
	vector<uint64_t> masks;

public:
	FragmentFactory();
	uint64_t coverage(int x1, int x2, int y1, int y2, int z1, int z2) const;
};

FragmentFactory::FragmentFactory()
{
	masks.resize(48);

	for (int k = 0; k < 16; k++)
	{
		int b = k % 4;
		int a = k / 4;

		uint64_t one = 1;
		uint64_t x_segment = 0;
		uint64_t y_segment = 0;
		uint64_t z_segment = 0;

		for (int l = a; l <= b; l++)
			for (int m = 0; m < 16; m++)
			{
				x_segment |= one << ((l / 2 % 2 + m / 8 % 2 * 2) * 16 + (l % 2 + m / 4 % 2 * 2) * 4 + (m % 4));
				y_segment |= one << ((l / 2 % 2 * 2 + m / 8 % 2) * 16 + (l % 2 * 2 + m / 4 % 2) * 4 + (m % 4));
				z_segment |= one << (m * 4 + l);
			}

		masks[0 + k] = x_segment;
		masks[16 + k] = y_segment;
		masks[32 + k] = z_segment;
	}
}

uint64_t FragmentFactory::coverage(int x1, int x2, int y1, int y2, int z1, int z2) const
{
	uint64_t x_segment = masks[0 + x1 * 4 + x2];
	uint64_t y_segment = masks[16 + y1 * 4 + y2];
	uint64_t z_segment = masks[32 + z1 * 4 + z2];
	uint64_t coverage = x_segment & y_segment & z_segment;
	return coverage;
}

template<typename T>
T clamp(T v, T lb, T ub)
{
	return std::min(std::max(v, lb), ub);
}

// Bounds computation utilities, similar to GPUQuad.hlsl
void UpdateClipRegionRoot(float nc,          // Tangent plane x/y normal coordinate (view space)
	float lc,          // Light x/y coordinate (view space)
	float lz,          // Light z coordinate (view space)
	float lightRadius,
	float cameraScale, // Project scale for coordinate (_11 or _22 for x/y respectively)
	float& clipMin,
	float& clipMax)
{
	float nz = (lightRadius - nc * lc) / lz;
	float pz = (lc * lc + lz * lz - lightRadius * lightRadius) / (lz - (nz / nc) * lc);

	if (pz > 0.0f) {
		float c = -nz * cameraScale / nc;
		if (nc > 0.0f)
		{                      // Left side boundary
			clipMin = std::max(clipMin, c);
		}
		else
		{                       // Right side boundary
			clipMax = std::min(clipMax, c);
		}
	}
}

void UpdateClipRegion(float lc,          // Light x/y coordinate (view space)
	float lz,          // Light z coordinate (view space)
	float lightRadius,
	float cameraScale, // Project scale for coordinate (_11 or _22 for x/y respectively)
	float& clipMin,
	float& clipMax)
{
	float rSq = lightRadius * lightRadius;
	float lcSqPluslzSq = lc * lc + lz * lz;
	float d = rSq * lc * lc - lcSqPluslzSq * (rSq - lz * lz);

	if (d > 0)
	{
		float a = lightRadius * lc;
		float b = sqrtf(d);
		float nx0 = (a + b) / lcSqPluslzSq;
		float nx1 = (a - b) / lcSqPluslzSq;

		UpdateClipRegionRoot(nx0, lc, lz, lightRadius, cameraScale, clipMin, clipMax);
		UpdateClipRegionRoot(nx1, lc, lz, lightRadius, cameraScale, clipMin, clipMax);
	}
}

// Returns bounding box [min.xy, max.xy] in clip [-1, 1] space.
DirectX::XMFLOAT4 ComputeClipRegion(DirectX::XMFLOAT3 lightPosView, float lightRadius,
	const DirectX::XMFLOAT4X4& mCameraProj, const DirectX::XMFLOAT4& mCameraNearFar)
{
	// Early out with empty rectangle if the light is too far behind the view frustum
	DirectX::XMFLOAT4 clipRegion = DirectX::XMFLOAT4(1, 1, 0, 0);
	if (lightPosView.z + lightRadius >= mCameraNearFar.x) {
		DirectX::XMFLOAT2 clipMin = DirectX::XMFLOAT2(-1.0f, -1.0f);
		DirectX::XMFLOAT2 clipMax = DirectX::XMFLOAT2(1.0f, 1.0f);

		UpdateClipRegion(lightPosView.x, lightPosView.z, lightRadius, mCameraProj._11, clipMin.x, clipMax.x);
		UpdateClipRegion(-lightPosView.y, lightPosView.z, lightRadius, mCameraProj._22, clipMin.y, clipMax.y);

		clipRegion = DirectX::XMFLOAT4(clipMin.x, clipMin.y, clipMax.x, clipMax.y);
	}

	return clipRegion;
}

void GenerateLightFragments(const FragmentFactory& fragmentFactory, LightGridBuilder* builder,
	const CFirstPersonCamera* viewerCamera, PointLight* light, int lightIndex)
{
	LightGridDimensions dim = builder->dimensions();
	DirectX::XMFLOAT4 mCameraNearFar = DirectX::XMFLOAT4(viewerCamera->GetFarClip(), viewerCamera->GetNearClip(), 0.0f, 0.0f);

	DirectX::XMFLOAT4X4 mCameraProj;
	DirectX::XMStoreFloat4x4(&mCameraProj, viewerCamera->GetProjMatrix());

	// compute view space quad
	DirectX::XMFLOAT4 clipRegion = ComputeClipRegion(light->positionView, light->attenuationEnd, mCameraProj, mCameraNearFar);
	// map coordinates to [0..1]
	clipRegion = DirectX::XMFLOAT4(0.5f * (clipRegion.x + 1.0f), 0.5f * (clipRegion.y + 1.0f), 0.5f * (clipRegion.z + 1.0f), 0.5f * (clipRegion.w + 1.0f));

	int intClipRegion[4];
	intClipRegion[0] = (int)(clipRegion.x * dim.width);
	intClipRegion[1] = (int)(clipRegion.y * dim.height);
	intClipRegion[2] = (int)(clipRegion.z * dim.width);
	intClipRegion[3] = (int)(clipRegion.w * dim.height);

	if (intClipRegion[0] < 0) intClipRegion[0] = 0;
	if (intClipRegion[1] < 0) intClipRegion[1] = 0;
	if (intClipRegion[2] >= dim.width) intClipRegion[2] = dim.width - 1;
	if (intClipRegion[3] >= dim.height) intClipRegion[3] = dim.height - 1;

	float center_z = (light->positionView.z - mCameraNearFar.x) / (mCameraNearFar.y - mCameraNearFar.x);
	float dist_z = light->attenuationEnd / (mCameraNearFar.y - mCameraNearFar.x);

	int intZBounds[2];
	intZBounds[0] = (int)((center_z - dist_z) * dim.depth);
	intZBounds[1] = (int)((center_z + dist_z) * dim.depth);

	if (intZBounds[0] < 0) intZBounds[0] = 0;
	if (intZBounds[1] >= dim.depth) intZBounds[1] = dim.depth - 1;


	for (int y = intClipRegion[1] / 4; y <= intClipRegion[3] / 4; y++)
		for (int x = intClipRegion[0] / 4; x <= intClipRegion[2] / 4; x++)
			for (int z = intZBounds[0] / 4; z <= intZBounds[1] / 4; z++)
			{
				int x1 = clamp(intClipRegion[0] - x * 4, 0, 3);
				int x2 = clamp(intClipRegion[2] - x * 4, 0, 3);
				int y1 = clamp(intClipRegion[1] - y * 4, 0, 3);
				int y2 = clamp(intClipRegion[3] - y * 4, 0, 3);
				int z1 = clamp(intZBounds[0] - z * 4, 0, 3);
				int z2 = clamp(intZBounds[1] - z * 4, 0, 3);

				uint64_t coverage = 0;
				coverage = fragmentFactory.coverage(x1, x2, y1, y2, z1, z2);

				builder->pushFragment(dim.cellIndex(x, y, z), lightIndex, coverage);
			}
}

void RasterizeLights(LightGridBuilder* builder, const CFirstPersonCamera* viewerCamera, PointLight lights[], int lightCount)
{
	ispc::Camera camera;
	FragmentFactory fragmentFactory;
	DirectX::XMFLOAT4X4 mCameraProj;
	DirectX::XMStoreFloat4x4(&mCameraProj, viewerCamera->GetProjMatrix());

	// z is flipped...
	camera.m_far = viewerCamera->GetNearClip();
	camera.m_near = viewerCamera->GetFarClip();
	camera.proj11 = mCameraProj._11;
	camera.proj22 = mCameraProj._22;

	bool simd = true;

	if (simd)
	{
		vector<ispc::LightBounds> bounds;
		bounds.resize(lightCount);
		LightGridDimensions dim = builder->dimensions();
		ispc::LightGridDimensions* pdim = (ispc::LightGridDimensions*)&dim;

		ispc::CoarseRasterizeLights((ispc::PointLight*)lights, &bounds[0], lightCount, &camera, pdim);

		// static to avoid large re-allocations
		static vector<ispc::Fragment> fragments;
		fragments.resize(0);

		for (int lightIndex = 0; lightIndex < lightCount; lightIndex++)
		{
			ispc::LightBounds region = bounds[lightIndex];

			for (int y = region.p1[1] / 4; y <= region.p2[1] / 4; y++)
				for (int x = region.p1[0] / 4; x <= region.p2[0] / 4; x++)
					for (int z = region.p1[2] / 4; z <= region.p2[2] / 4; z++)
					{
						ispc::Fragment fragment;
						fragment.cellIndex = dim.cellIndex(x, y, z);
						fragment.lightIndex = lightIndex;
						fragments.push_back(fragment);
					}
		}

		int fragCount = (int)fragments.size();
		ispc::FineRasterizeLights((ispc::PointLight*)lights, &fragments[0], fragCount, &camera, pdim);

		for (int fragIndex = 0; fragIndex < fragCount; fragIndex++)
		{
			ispc::Fragment fragment = fragments[fragIndex];

			builder->pushFragment(fragment.cellIndex, fragment.lightIndex, fragment.coverage);
		}
	}
	else
	{
		// warning: scalar version does coarser (AABB) culling
		for (int lightIndex = 0; lightIndex < lightCount; lightIndex++)
		{
			GenerateLightFragments(fragmentFactory, builder, viewerCamera, &lights[lightIndex], lightIndex);
		}
	}
}
