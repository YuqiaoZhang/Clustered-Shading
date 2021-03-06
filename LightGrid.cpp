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

#include <algorithm>
#include <cassert>
#include "LightGrid.h"

LightGridBuilder::LightGridBuilder() : dim(0, 0, 0)
{
	// nothing
}

void LightGridBuilder::reset(LightGridDimensions _dim)
{
	dim = _dim;
	int cellCount = dim.width * dim.height * dim.depth / 64;
	lightIndexLists.resize(cellCount);
	coverageLists.resize(cellCount);
}

LightGridDimensions LightGridBuilder::dimensions()
{
	return dim;
}

size_t LightGridBuilder::cellCount()
{
	assert(lightIndexLists.size() == coverageLists.size());
	return lightIndexLists.size();
}

void LightGridBuilder::clearAllFragments()
{
	for (size_t cellIndex = 0; cellIndex < cellCount(); cellIndex++)
	{
		lightIndexLists[cellIndex].resize(0);
		coverageLists[cellIndex].resize(0);
	}
}

void LightGridBuilder::pushFragment(int cellIndex, int lightIndex, uint64_t coverage)
{
	if (coverage == 0) return;
	lightIndexLists[cellIndex].push_back(lightIndex);
	coverageLists[cellIndex].push_back(coverage);
}

void LightGridBuilder::buildAndUpload(void* gpuBuffer, size_t bufferSize)
{
	tempBuffer.resize(bufferSize + 256 * 1024); // 256 KB margin: max allocation per cell

	size_t headerBytes = cellCount() * 64 * 16; // uint4: 16 bytes per entry
	allocatedBytes = headerBytes;
	size_t uploadedBytes = headerBytes;

	for (int y = 0; y < dim.height / 4; y++)
		for (int x = 0; x < dim.width / 4; x++)
			for (int z = 0; z < dim.depth / 4; z++)
			{
				buildFlatEntries(x, y, z);

				if (gpuBuffer) memcpy(((char*)gpuBuffer) + uploadedBytes, &tempBuffer[uploadedBytes], allocatedBytes - uploadedBytes);
				uploadedBytes = allocatedBytes;

				assert(allocatedBytes <= bufferSize && "gpu buffer not big enough");
			}

	if (gpuBuffer) memcpy(gpuBuffer, &tempBuffer[0], headerBytes);
}

inline int getFineIndex(int xx, int yy)
{
	static const int fineIndexTable[][4] =
	{
		{ 0, 1, 4, 5 },
		{ 2, 3, 6, 7 },
		{ 8, 9, 12, 13 },
		{ 10, 11, 14, 15 },
	};
	return fineIndexTable[yy][xx];
}

uint32_t swap_word_pair(uint32_t pair)
{
	return (pair << 16) | (pair >> 16);
}

void LightGridBuilder::buildFlatEntries(int _x, int _y, int _z)
{
	int cellIndex = dim.cellIndex(_x, _y, _z);
	vector<int>& lightIndexList = lightIndexLists[cellIndex];
	vector<uint64_t>& coverageList = coverageLists[cellIndex];
	size_t count = lightIndexList.size();
	assert(count == coverageList.size());

	if (count == 0)
	{
		for (int entryIndex = 0; entryIndex < 64; entryIndex++)
		{
			int yy = entryIndex / 16;
			int xx = entryIndex / 4 % 4;
			int zz = entryIndex % 4;

			int x = _x * 4 + xx;
			int y = _y * 4 + yy;
			int z = _z * 4 + zz;

			int headerIndex = (y * dim.width + x) * dim.depth + z;
			uint32_t* entry_ptr = (uint32_t*)&tempBuffer[16 * headerIndex];

			entry_ptr[0] = 0; // list size: 0
		}

		return;
	}

	int* lightIndexList_ptr = &lightIndexList[0];
	uint32_t* coverageList_ptr = (uint32_t*)&coverageList[0];

	for (int entryIndex = 0; entryIndex < 64; entryIndex++)
	{
		int yy = entryIndex / 16;
		int xx = entryIndex / 4 % 4;
		int zz = entryIndex % 4;

		int x = _x * 4 + xx;
		int y = _y * 4 + yy;
		int z = _z * 4 + zz;

		size_t headerIndex = (y * dim.width + x) * dim.depth + z;
		uint32_t* entry_ptr = (uint32_t*)&tempBuffer[16 * headerIndex];
		uint16_t* tail_ptr = (uint16_t*)&tempBuffer[allocatedBytes];

		int fineIndex = getFineIndex(xx, yy) * 4 + zz;
		uint64_t mask = uint64_t(1) << fineIndex;
		uint32_t sub_mask = (uint32_t)mask;
		uint32_t* sub_coverageList_ptr = coverageList_ptr;
		if (sub_mask == 0)
		{
			sub_mask = mask >> 32;
			sub_coverageList_ptr++;
		}

		size_t cursor = 0;
		for (size_t k = 0; k < count; k++)
		{
			tail_ptr[cursor] = uint16_t(lightIndexList_ptr[k]);
			cursor += !!(sub_coverageList_ptr[k * 2] & sub_mask);
		}

		entry_ptr[1] = swap_word_pair(*(uint32_t*)&tail_ptr[cursor - 2]);
		entry_ptr[2] = swap_word_pair(*(uint32_t*)&tail_ptr[cursor - 4]);
		entry_ptr[3] = swap_word_pair(*(uint32_t*)&tail_ptr[cursor - 6]);

		size_t list_size = cursor;
		assert(list_size < 0x100);
		assert(allocatedBytes / 16 < 0x1000000);
		entry_ptr[0] = uint32_t(((allocatedBytes / 16) << 8) | list_size);

		allocatedBytes += (std::max(size_t(6), list_size) - 6 + 7) / 8 * 16;
	}
}
