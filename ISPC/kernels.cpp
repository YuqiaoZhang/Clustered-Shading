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

#define WIN32_LEAN_AND_MEAN 1
#define NOMINMAX 1
#include <Windows.h>
#include <DirectXMath.h>
#include <stdint.h>
#include <algorithm>

#include "kernels.h"

namespace ispc
{
	// Bounds computation utilities, similar to FragmentFactory.cpp
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
			float b = sqrt(d);
			float nx0 = (a + b) / lcSqPluslzSq;
			float nx1 = (a - b) / lcSqPluslzSq;

			UpdateClipRegionRoot(nx0, lc, lz, lightRadius, cameraScale, clipMin, clipMax);
			UpdateClipRegionRoot(nx1, lc, lz, lightRadius, cameraScale, clipMin, clipMax);
		}
	}

	// Returns bounding box [min.xy, max.xy] in clip [-1, 1] space.
	DirectX::XMFLOAT4 ComputeClipRegion(DirectX::XMFLOAT3 lightPosView, float lightRadius, Camera camera[])
	{
		// Early out with empty rectangle if the light is too far behind the view frustum
		DirectX::XMFLOAT4 clipRegion = { 1, 1, 0, 0 };
		if (lightPosView.z + lightRadius >= camera->m_near) {
			DirectX::XMFLOAT2 clipMin = { -1.0f, -1.0f };
			DirectX::XMFLOAT2 clipMax = { 1.0f, 1.0f };

			UpdateClipRegion(lightPosView.x, lightPosView.z, lightRadius, camera->proj11, clipMin.x, clipMax.x);
			UpdateClipRegion(-lightPosView.y, lightPosView.z, lightRadius, camera->proj22, clipMin.y, clipMax.y);

			DirectX::XMFLOAT4 v = { clipMin.x, clipMin.y, clipMax.x, clipMax.y };
			clipRegion = v;
		}

		return clipRegion;
	}

	DirectX::XMFLOAT3 to_float3(float v[3])
	{
		DirectX::XMFLOAT3 r = { v[0], v[1], v[2] };
		return r;
	}

	void GenerateLightBounds(PointLight light, LightBounds& box, Camera camera[], LightGridDimensions dim[])
	{
		// compute view space quad
		DirectX::XMFLOAT4 clipRegion = ComputeClipRegion(to_float3(light.positionView), light.attenuationEnd, camera);

		// map coordinates to [0..1]
		clipRegion.x = (clipRegion.x + 1.0f) / 2;
		clipRegion.y = (clipRegion.y + 1.0f) / 2;
		clipRegion.z = (clipRegion.z + 1.0f) / 2;
		clipRegion.w = (clipRegion.w + 1.0f) / 2;

		int intClipRegion[4];
		intClipRegion[0] = (int)(clipRegion.x * dim->width);
		intClipRegion[1] = (int)(clipRegion.y * dim->height);
		intClipRegion[2] = (int)(clipRegion.z * dim->width);
		intClipRegion[3] = (int)(clipRegion.w * dim->height);

		if (intClipRegion[0] < 0) intClipRegion[0] = 0;
		if (intClipRegion[1] < 0) intClipRegion[1] = 0;
		if (intClipRegion[2] >= dim->width) intClipRegion[2] = dim->width - 1;
		if (intClipRegion[3] >= dim->height) intClipRegion[3] = dim->height - 1;

		float center_z = (light.positionView[2] - camera->m_near) / (camera->m_far - camera->m_near);
		float dist_z = light.attenuationEnd / (camera->m_far - camera->m_near);

		int intZBounds[2];
		intZBounds[0] = (int)((center_z - dist_z) * dim->depth);
		intZBounds[1] = (int)((center_z + dist_z) * dim->depth);

		if (intZBounds[0] < 0) intZBounds[0] = 0;
		if (intZBounds[1] >= dim->depth) intZBounds[1] = dim->depth - 1;

		box.p1[0] = intClipRegion[0];
		box.p2[0] = intClipRegion[2];
		box.p1[1] = intClipRegion[1];
		box.p2[1] = intClipRegion[3];

		box.p1[2] = intZBounds[0];
		box.p2[2] = intZBounds[1];
	}

	void CoarseRasterizeLights(PointLight lights[], LightBounds bounds[], int lightCount, Camera camera[], LightGridDimensions dim[])
	{
		for (int idx = 0; idx < lightCount; ++idx)
		{
			LightBounds box;
			GenerateLightBounds(lights[idx], box, camera, dim);
			bounds[idx] = box;
		}
	}

	float dot(DirectX::XMFLOAT3 a, DirectX::XMFLOAT3 b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	float length(DirectX::XMFLOAT3 v)
	{
		return sqrt(dot(v, v));
	}

	uint64_t ComputeCoverage(int cellIndex, DirectX::XMFLOAT3 lightPosition, float lightSize, Camera camera[], LightGridDimensions dim[])
	{
		int cz = cellIndex % (dim->depth / 4);
		int cx = (cellIndex / (dim->depth / 4)) % (dim->width / 4);
		int cy = cellIndex / (dim->depth / 4 * dim->width / 4);

		uint64_t coverage = 0;
		for (int zz = 0; zz <= 3; zz++)
			for (int yy = 0; yy <= 3; yy++)
				for (int xx = 0; xx <= 3; xx++)
				{
					int fineIndex = (yy / 2 % 2) * 32 + (xx / 2 % 2) * 16 + (yy % 2) * 8 + (xx % 2) * 4 + (zz % 4);

					int x = cx * 4 + xx;
					int y = cy * 4 + yy;
					int z = cz * 4 + zz;

					// sub-frustrum bounds in view space        
					float minZ = (z - 0) * 1.0f / dim->depth * (camera->m_far - camera->m_near) + camera->m_near;
					float maxZ = (z + 1) * 1.0f / dim->depth * (camera->m_far - camera->m_near) + camera->m_near;

					float minZminX = -(1 - 2.0f / dim->width * (x - 0)) * minZ / camera->proj11;
					float minZmaxX = -(1 - 2.0f / dim->width * (x + 1)) * minZ / camera->proj11;
					float minZminY = (1 - 2.0f / dim->height * (y - 0)) * minZ / camera->proj22;
					float minZmaxY = (1 - 2.0f / dim->height * (y + 1)) * minZ / camera->proj22;

					float maxZminX = -(1 - 2.0f / dim->width * (x - 0)) * maxZ / camera->proj11;
					float maxZmaxX = -(1 - 2.0f / dim->width * (x + 1)) * maxZ / camera->proj11;
					float maxZminY = (1 - 2.0f / dim->height * (y - 0)) * maxZ / camera->proj22;
					float maxZmaxY = (1 - 2.0f / dim->height * (y + 1)) * maxZ / camera->proj22;

					// heuristic plane separation test - works pretty well in practice
					DirectX::XMFLOAT3 minZcenter = { (minZminX + minZmaxX) / 2, (minZminY + minZmaxY) / 2, minZ };
					DirectX::XMFLOAT3 maxZcenter = { (maxZminX + maxZmaxX) / 2, (maxZminY + maxZmaxY) / 2, maxZ };

					DirectX::XMFLOAT3 center;
					DirectX::XMStoreFloat3(&center, DirectX::XMVectorScale(DirectX::XMVectorAdd(DirectX::XMLoadFloat3(&minZcenter), DirectX::XMLoadFloat3(&maxZcenter)), 0.5));

					DirectX::XMFLOAT3 normal;
					DirectX::XMStoreFloat3(&normal, DirectX::XMVector3Normalize(DirectX::XMVectorSubtract(DirectX::XMLoadFloat3(&center), DirectX::XMLoadFloat3(&lightPosition))));

					// compute minimum distance of all 8 corners to the tangent plane, with a few shortcuts (saves 14 muls)
					float min_d1 = -dot(normal, lightPosition);
					float min_d2 = min_d1;
					min_d1 += std::min(normal.x * minZminX, normal.x * minZmaxX);
					min_d1 += std::min(normal.y * minZminY, normal.y * minZmaxY);
					min_d1 += normal.z * minZ;
					min_d2 += std::min(normal.x * maxZminX, normal.x * maxZmaxX);
					min_d2 += std::min(normal.y * maxZminY, normal.y * maxZmaxY);
					min_d2 += normal.z * maxZ;
					float min_d = std::min(min_d1, min_d2);
					bool separated = min_d > lightSize;

					uint64_t one = 1;
					if (!separated)
						coverage |= one << fineIndex;
				}

		return coverage;
	}

	float sq(float v)
	{
		return v * v;
	}

	inline  int getFineIndex(int xx, int yy)
	{
		static const  int fineIndexTable[][4] =
		{
			{ 0, 1, 4, 5 },
			{ 2, 3, 6, 7 },
			{ 8, 9, 12, 13 },
			{ 10, 11, 14, 15 },
		};
		return fineIndexTable[yy][xx];
	}

	// same as ComputeCoverage, but faster/uglier
	uint64_t ComputeCoverage2(int cellIndex, DirectX::XMFLOAT3 lightPosition, float lightSize, Camera camera[], LightGridDimensions dim[])
	{
		int cz = cellIndex % (dim->depth / 4);
		int cx = (cellIndex / (dim->depth / 4)) % (dim->width / 4);
		int cy = cellIndex / (dim->depth / 4 * dim->width / 4);

		uint64_t coverage = 0;
		for (int zz = 0; zz <= 3; zz++)
		{
			// Z
			int z = cz * 4 + zz;
			float minZ = (z - 0) * 1.0f / dim->depth * (camera->m_far - camera->m_near) + camera->m_near;
			float maxZ = (z + 1) * 1.0f / dim->depth * (camera->m_far - camera->m_near) + camera->m_near;

			float centerZ = (minZ + maxZ) * 0.5;
			float normalZ = centerZ - lightPosition.z;

			float d0Z = normalZ * lightPosition.z;
			float min_d1Z = -d0Z + normalZ * minZ;
			float min_d2Z = -d0Z + normalZ * maxZ;

			// X
			float minZmulX = 2.0f / dim->width * minZ / camera->proj11;
			float minZaddX = -minZ / camera->proj11;
			float maxZmulX = 2.0f / dim->width * maxZ / camera->proj11;
			float maxZaddX = -maxZ / camera->proj11;

			float min_d1X[4];
			float min_d2X[4];
			float normal2X[4];

			for (int xx = 0; xx <= 3; xx++)
			{
				int x = cx * 4 + xx;

				float minZminX = (x - 0) * minZmulX + minZaddX;
				float minZmaxX = (x + 1) * minZmulX + minZaddX;
				float maxZminX = (x - 0) * maxZmulX + maxZaddX;
				float maxZmaxX = (x + 1) * maxZmulX + maxZaddX;

				float centerX = (minZminX + minZmaxX + maxZminX + maxZmaxX) * 0.25;
				float normalX = centerX - lightPosition.x;

				float d0X = normalX * lightPosition.x;
				min_d1X[xx] = -d0X + std::min(normalX * minZminX, normalX * minZmaxX);
				min_d2X[xx] = -d0X + std::min(normalX * maxZminX, normalX * maxZmaxX);
				normal2X[xx] = sq(normalX);
			}

			// Y
			float minZmulY = -2.0f / dim->height * minZ / camera->proj22;
			float minZaddY = minZ / camera->proj22;
			float maxZmulY = -2.0f / dim->height * maxZ / camera->proj22;
			float maxZaddY = maxZ / camera->proj22;

			float min_d1Y[4];
			float min_d2Y[4];
			float normal2Y[4];

			for (int yy = 0; yy <= 3; yy++)
			{
				int y = cy * 4 + yy;

				float minZminY = (y - 0) * minZmulY + minZaddY;
				float minZmaxY = (y + 1) * minZmulY + minZaddY;
				float maxZminY = (y - 0) * maxZmulY + maxZaddY;
				float maxZmaxY = (y + 1) * maxZmulY + maxZaddY;

				float centerY = (minZminY + minZmaxY + maxZminY + maxZmaxY) * 0.25;
				float normalY = centerY - lightPosition.y;

				float d0Y = normalY * lightPosition.y;
				min_d1Y[yy] = -d0Y + std::min(normalY * minZminY, normalY * minZmaxY);
				min_d2Y[yy] = -d0Y + std::min(normalY * maxZminY, normalY * maxZmaxY);
				normal2Y[yy] = sq(normalY);
			}

			// rasterize a Z slice
			for (int yy = 0; yy <= 3; yy++)
			{
				for (int xx = 0; xx <= 3; xx++)
				{
					unsigned int fineIndex = getFineIndex(xx, yy) * 4 + zz;

					float normal2 = sq(normalZ);
					normal2 += normal2X[xx];
					normal2 += normal2Y[yy];

					float min_d1 = min_d1X[xx] + min_d1Y[yy] + min_d1Z;
					float min_d2 = min_d2X[xx] + min_d2Y[yy] + min_d2Z;
					float min_d = std::min(min_d1, min_d2);

					bool separated = sq(min_d) > sq(lightSize) * normal2;
					if (min_d < 0) separated = false;

					uint64_t one = 1;
					if (!separated)
						coverage |= one << fineIndex;
				}
			}
		}

		return coverage;
	}

	void FineRasterizeLights(PointLight lights[], Fragment fragments[], int fragmentCount, Camera camera[], LightGridDimensions dim[])
	{
		for (int idx = 0; idx < fragmentCount; ++idx)
		{
			PointLight light = lights[fragments[idx].lightIndex];
			fragments[idx].coverage = ComputeCoverage2(fragments[idx].cellIndex, to_float3(light.positionView), light.attenuationEnd, camera, dim);
			//fragments[idx].coverage = 0xFFFFFFFFFFFFFFFF;
		}
	}
}
