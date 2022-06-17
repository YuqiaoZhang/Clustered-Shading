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

#ifndef _CLUSTERED_HLSLI_
#define _CLUSTERED_HLSLI_ 1

#include "Rendering.hlsli"

StructuredBuffer<uint4> gLightGrid : register(t6);

void fill_array4(out int array[4], int4 src)
{
	array[0] = src.x;
	array[1] = src.y;
	array[2] = src.z;
	array[3] = src.w;
}

float4 ClusteredPS(GeometryVSOut input) : SV_Target
{
	uint3 grid = uint3(2, 1, 8) * mUI.clusteredGridScale;

	float3 lit = float3(0.0f, 0.0f, 0.0f);
	SurfaceData surface = ComputeSurfaceDataFromGeometry(input);

	float2 screenPosition = input.position.xy / mFramebufferDimensions.xy;
	float zPosition = (surface.positionView.z - mCameraNearFar.x) / (mCameraNearFar.y - mCameraNearFar.x);

	uint3 clusterPosition = uint3(screenPosition * grid.xy, zPosition * grid.z);
	uint cluster_index = (clusterPosition.y * grid.x + clusterPosition.x) * grid.z + clusterPosition.z;

	int lightIndexBlock[4];
	fill_array4(lightIndexBlock, gLightGrid[cluster_index]);

	int list_size = lightIndexBlock[0] & 255;
	int list_index = lightIndexBlock[0] >> 8;
	int light_count = list_size;

	[loop]
	for (int k = 2; k < list_size + 2; k++)
	{
		int lightIndex = (lightIndexBlock[(k & 7) >> 1] >> ((k & 1) << 4)) & 0xFFFF;
		if ((k & 7) == 7) fill_array4(lightIndexBlock, gLightGrid[list_index++]);

		PointLight light = gLight[lightIndex];
		AccumulateBRDF(surface, light, lit);
	}

	[flatten] if (mUI.visualizeLightCount)
	{
		lit = (float(light_count) * rcp(255.0f)).xxx;
	}

	return float4(lit, 1.0f);
}

#endif // CLUSTERED_HLSL
