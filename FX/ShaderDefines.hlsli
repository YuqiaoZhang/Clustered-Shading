// Copyright 2010 Intel Corporation
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

#ifndef _SHADERDEFINES_HLSLI_
#define _SHADERDEFINES_HLSLI_ 1

#define MAX_LIGHTS_POWER 12
#define MAX_LIGHTS (1<<MAX_LIGHTS_POWER)

// reduce maximum light list size per tile for better occupancy (more realistic performance)
#define MAX_SMEM_LIGHTS 512

// This determines the tile size for light binning and associated tradeoffs
#define COMPUTE_SHADER_TILE_GROUP_DIM 16
#define COMPUTE_SHADER_TILE_GROUP_SIZE (COMPUTE_SHADER_TILE_GROUP_DIM*COMPUTE_SHADER_TILE_GROUP_DIM)

// If enabled, defers scheduling of per-sample-shaded pixels until after sample 0
// has been shaded across the whole tile. This allows better SIMD packing and scheduling.
// This should basically always be left enabled in practice since it's faster everywhere
// but we maintain the legacy path for benckmarking comparison for now.
#define DEFER_PER_SAMPLE 1

// If enabled, add a plane test heuristic to the frustum test
#define MORE_CULLING 0

#endif
