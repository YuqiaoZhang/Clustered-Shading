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

#include "App.h"
#include "ColorUtil.h"
#include "FX/ShaderDefines.hlsli"
#include "LightGrid.h"
#include "FragmentFactory.h"
#include <limits>
#include <sstream>
#include <random>
#include <algorithm>

#include "FX/Rendering_GeometryVS.hlsl.h"
#include "FX/GBuffer_GBufferPS.hlsl.h"
#include "FX/GBuffer_GBufferAlphaTestPS.hlsl.h"
#include "FX/Forward_ForwardPS.hlsl.h"
#include "FX/Forward_ForwardAlphaTestPS.hlsl.h"
#include "FX/Forward_ForwardAlphaTestOnlyPS.hlsl.h"
#include "FX/Clustered_ClusteredPS.hlsl.h"
#include "FX/Rendering_FullScreenTriangleVS.hlsl.h"
#include "FX/SkyboxToneMap_SkyboxVS.hlsl.h"
#include "FX/SkyboxToneMap_SkyboxPS.hlsl.h"
#include "FX/GBuffer_RequiresPerSampleShadingPS.hlsl.h"
#include "FX/BasicLoop_BasicLoopPS.hlsl.h"
#include "FX/BasicLoop_BasicLoopPerSamplePS.hlsl.h"
#include "FX/ComputeShaderTile_ComputeShaderTileCS.hlsl.h"
#include "FX/GPUQuad_GPUQuadVS.hlsl.h"
#include "FX/GPUQuad_GPUQuadGS.hlsl.h"
#include "FX/GPUQuad_GPUQuadPS.hlsl.h"
#include "FX/GPUQuad_GPUQuadPerSamplePS.hlsl.h"
#include "FX/GPUQuadDL_GPUQuadDLPS.hlsl.h"
#include "FX/GPUQuadDL_GPUQuadDLPerSamplePS.hlsl.h"
#include "FX/GPUQuadDL_GPUQuadDLResolvePS.hlsl.h"
#include "FX/GPUQuadDL_GPUQuadDLResolvePerSamplePS.hlsl.h"

using std::shared_ptr;

// NOTE: Must match layout of shader constant buffers
__declspec(align(16)) struct PerFrameConstants
{
	DirectX::XMFLOAT4X4 mCameraWorldViewProj;
	DirectX::XMFLOAT4X4 mCameraWorldView;
	DirectX::XMFLOAT4X4 mCameraViewProj;
	DirectX::XMFLOAT4X4 mCameraProj;
	DirectX::XMFLOAT4 mCameraNearFar;

	unsigned int mFramebufferDimensionsX;
	unsigned int mFramebufferDimensionsY;
	unsigned int mFramebufferDimensionsZ;
	unsigned int mFramebufferDimensionsW;

	UIConstants mUI;
};

App::App(ID3D11Device* d3dDevice, double log2ActiveLights, unsigned int msaaSamples)
	: mMSAASamples(msaaSamples), mTotalTime(0.0f), mActiveLights(0), mLightBuffer(0), mLightGridBuffer(0), mDepthBufferReadOnlyDSV(0)
{
	HRESULT hr;

	std::string msaaSamplesStr;
	{
		std::ostringstream oss;
		oss << mMSAASamples;
		msaaSamplesStr = oss.str();
	}

	// Create shaders
	hr = d3dDevice->CreateVertexShader(Rendering_GeometryVS, sizeof(Rendering_GeometryVS), NULL, &mGeometryVS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(GBuffer_GBufferPS, sizeof(GBuffer_GBufferPS), NULL, &mGBufferPS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(GBuffer_GBufferAlphaTestPS, sizeof(GBuffer_GBufferAlphaTestPS), NULL, &mGBufferAlphaTestPS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(Forward_ForwardPS, sizeof(Forward_ForwardPS), NULL, &mForwardPS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(Forward_ForwardAlphaTestPS, sizeof(Forward_ForwardAlphaTestPS), NULL, &mForwardAlphaTestPS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(Forward_ForwardAlphaTestOnlyPS, sizeof(Forward_ForwardAlphaTestOnlyPS), NULL, &mForwardAlphaTestOnlyPS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(Clustered_ClusteredPS, sizeof(Clustered_ClusteredPS), NULL, &mClusteredPS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreateVertexShader(Rendering_FullScreenTriangleVS, sizeof(Rendering_FullScreenTriangleVS), NULL, &mFullScreenTriangleVS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreateVertexShader(SkyboxToneMap_SkyboxVS, sizeof(SkyboxToneMap_SkyboxVS), NULL, &mSkyboxVS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(SkyboxToneMap_SkyboxPS, sizeof(SkyboxToneMap_SkyboxPS), NULL, &mSkyboxPS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(GBuffer_RequiresPerSampleShadingPS, sizeof(GBuffer_RequiresPerSampleShadingPS), NULL, &mRequiresPerSampleShadingPS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(BasicLoop_BasicLoopPS, sizeof(BasicLoop_BasicLoopPS), NULL, &mBasicLoopPS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(BasicLoop_BasicLoopPS, sizeof(BasicLoop_BasicLoopPS), NULL, &mBasicLoopPerSamplePS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreateComputeShader(ComputeShaderTile_ComputeShaderTileCS, sizeof(ComputeShaderTile_ComputeShaderTileCS), NULL, &mComputeShaderTileCS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreateVertexShader(GPUQuad_GPUQuadVS, sizeof(GPUQuad_GPUQuadVS), NULL, &mGPUQuadVS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreateGeometryShader(GPUQuad_GPUQuadGS, sizeof(GPUQuad_GPUQuadGS), NULL, &mGPUQuadGS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(GPUQuad_GPUQuadPS, sizeof(GPUQuad_GPUQuadPS), NULL, &mGPUQuadPS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(GPUQuad_GPUQuadPerSamplePS, sizeof(GPUQuad_GPUQuadPerSamplePS), NULL, &mGPUQuadPerSamplePS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(GPUQuadDL_GPUQuadDLPS, sizeof(GPUQuadDL_GPUQuadDLPS), NULL, &mGPUQuadDLPS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(GPUQuadDL_GPUQuadDLPerSamplePS, sizeof(GPUQuadDL_GPUQuadDLPerSamplePS), NULL, &mGPUQuadDLPerSamplePS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(GPUQuadDL_GPUQuadDLResolvePS, sizeof(GPUQuadDL_GPUQuadDLResolvePS), NULL, &mGPUQuadDLResolvePS);
	assert(SUCCEEDED(hr));

	hr = d3dDevice->CreatePixelShader(GPUQuadDL_GPUQuadDLResolvePerSamplePS, sizeof(GPUQuadDL_GPUQuadDLResolvePerSamplePS), NULL, &mGPUQuadDLResolvePerSamplePS);
	assert(SUCCEEDED(hr));

	// Create input layout
	{
		const D3D11_INPUT_ELEMENT_DESC layout[] =
		{
			{"position", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0},
			{"normal", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0},
			{"texCoord", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 24, D3D11_INPUT_PER_VERTEX_DATA, 0},
		};

		d3dDevice->CreateInputLayout(layout, ARRAYSIZE(layout), Rendering_GeometryVS, sizeof(Rendering_GeometryVS), &mMeshVertexLayout);
	}

	// Create standard rasterizer state
	{
		CD3D11_RASTERIZER_DESC desc(D3D11_DEFAULT);
		hr = d3dDevice->CreateRasterizerState(&desc, &mRasterizerState);
		assert(SUCCEEDED(hr));

		desc.CullMode = D3D11_CULL_NONE;
		hr = d3dDevice->CreateRasterizerState(&desc, &mDoubleSidedRasterizerState);
		assert(SUCCEEDED(hr));
	}

	{
		CD3D11_DEPTH_STENCIL_DESC desc(D3D11_DEFAULT);
		// NOTE: Complementary Z => GREATER test
		desc.DepthFunc = D3D11_COMPARISON_GREATER_EQUAL;
		hr = d3dDevice->CreateDepthStencilState(&desc, &mDepthState);
		assert(SUCCEEDED(hr));
	}

	// Stencil states for MSAA
	{
		CD3D11_DEPTH_STENCIL_DESC desc(
			FALSE, D3D11_DEPTH_WRITE_MASK_ZERO, D3D11_COMPARISON_GREATER_EQUAL,									   // Depth
			TRUE, 0xFF, 0xFF,																					   // Stencil
			D3D11_STENCIL_OP_REPLACE, D3D11_STENCIL_OP_REPLACE, D3D11_STENCIL_OP_REPLACE, D3D11_COMPARISON_ALWAYS, // Front face stencil
			D3D11_STENCIL_OP_REPLACE, D3D11_STENCIL_OP_REPLACE, D3D11_STENCIL_OP_REPLACE, D3D11_COMPARISON_ALWAYS  // Back face stencil
		);
		hr = d3dDevice->CreateDepthStencilState(&desc, &mWriteStencilState);
		assert(SUCCEEDED(hr));
	}
	{
		CD3D11_DEPTH_STENCIL_DESC desc(
			TRUE, D3D11_DEPTH_WRITE_MASK_ZERO, D3D11_COMPARISON_GREATER_EQUAL,							 // Depth
			TRUE, 0xFF, 0xFF,																			 // Stencil
			D3D11_STENCIL_OP_KEEP, D3D11_STENCIL_OP_KEEP, D3D11_STENCIL_OP_KEEP, D3D11_COMPARISON_EQUAL, // Front face stencil
			D3D11_STENCIL_OP_KEEP, D3D11_STENCIL_OP_KEEP, D3D11_STENCIL_OP_KEEP, D3D11_COMPARISON_EQUAL	 // Back face stencil
		);
		hr = d3dDevice->CreateDepthStencilState(&desc, &mEqualStencilState);
		assert(SUCCEEDED(hr));
	}

	// Create geometry phase blend state
	{
		CD3D11_BLEND_DESC desc(D3D11_DEFAULT);
		hr = d3dDevice->CreateBlendState(&desc, &mGeometryBlendState);
		assert(SUCCEEDED(hr));
	}

	// Create lighting phase blend state
	{
		CD3D11_BLEND_DESC desc(D3D11_DEFAULT);
		// Additive blending
		desc.RenderTarget[0].BlendEnable = true;
		desc.RenderTarget[0].SrcBlend = D3D11_BLEND_ONE;
		desc.RenderTarget[0].DestBlend = D3D11_BLEND_ONE;
		desc.RenderTarget[0].BlendOp = D3D11_BLEND_OP_ADD;
		desc.RenderTarget[0].SrcBlendAlpha = D3D11_BLEND_ONE;
		desc.RenderTarget[0].DestBlendAlpha = D3D11_BLEND_ONE;
		desc.RenderTarget[0].BlendOpAlpha = D3D11_BLEND_OP_ADD;
		hr = d3dDevice->CreateBlendState(&desc, &mLightingBlendState);
		assert(SUCCEEDED(hr));
	}

	// Create constant buffers
	{
		CD3D11_BUFFER_DESC desc(
			sizeof(PerFrameConstants),
			D3D11_BIND_CONSTANT_BUFFER,
			D3D11_USAGE_DYNAMIC,
			D3D11_CPU_ACCESS_WRITE);

		hr = d3dDevice->CreateBuffer(&desc, 0, &mPerFrameConstants);
		assert(SUCCEEDED(hr));
	}

	// Create sampler state
	{
		CD3D11_SAMPLER_DESC desc(D3D11_DEFAULT);
		desc.Filter = D3D11_FILTER_ANISOTROPIC;
		desc.AddressU = D3D11_TEXTURE_ADDRESS_WRAP;
		desc.AddressV = D3D11_TEXTURE_ADDRESS_WRAP;
		desc.AddressW = D3D11_TEXTURE_ADDRESS_WRAP;
		desc.MaxAnisotropy = 16;
		hr = d3dDevice->CreateSamplerState(&desc, &mDiffuseSampler);
		assert(SUCCEEDED(hr));
	}

	// Create skybox mesh
	hr = mSkyboxMesh.Create(d3dDevice, L"Media\\Skybox\\Skybox.sdkmesh");
	assert(SUCCEEDED(hr));

	// Create light grid buffer
	mLightGridBufferSize = 64 * 1024 * 32; // 32 MB
	if (!mLightGridBuffer)
	{
		mLightGridBuffer = new StructuredBuffer<LightGridEntry>(d3dDevice, mLightGridBufferSize, D3D11_BIND_SHADER_RESOURCE, true);
	}

	InitializeLightParameters(d3dDevice);
	SetActiveLights(d3dDevice, log2ActiveLights);
}

App::~App()
{
	mSkyboxMesh.Destroy();
	SAFE_RELEASE(mDepthBufferReadOnlyDSV);
	delete mLightBuffer;
	delete mLightGridBuffer;
	SAFE_RELEASE(mDiffuseSampler);
	SAFE_RELEASE(mPerFrameConstants);
	SAFE_RELEASE(mLightingBlendState);
	SAFE_RELEASE(mGeometryBlendState);
	SAFE_RELEASE(mEqualStencilState);
	SAFE_RELEASE(mWriteStencilState);
	SAFE_RELEASE(mDepthState);
	SAFE_RELEASE(mDoubleSidedRasterizerState);
	SAFE_RELEASE(mRasterizerState);
	SAFE_RELEASE(mMeshVertexLayout);
	SAFE_RELEASE(mSkyboxPS);
	SAFE_RELEASE(mSkyboxVS);
	SAFE_RELEASE(mComputeShaderTileCS);
	SAFE_RELEASE(mGPUQuadDLResolvePerSamplePS);
	SAFE_RELEASE(mGPUQuadDLResolvePS);
	SAFE_RELEASE(mGPUQuadDLPerSamplePS);
	SAFE_RELEASE(mGPUQuadDLPS);
	SAFE_RELEASE(mGPUQuadPerSamplePS);
	SAFE_RELEASE(mGPUQuadPS);
	SAFE_RELEASE(mGPUQuadGS);
	SAFE_RELEASE(mGPUQuadVS);
	SAFE_RELEASE(mRequiresPerSampleShadingPS);
	SAFE_RELEASE(mBasicLoopPerSamplePS);
	SAFE_RELEASE(mBasicLoopPS);
	SAFE_RELEASE(mFullScreenTriangleVS);
	SAFE_RELEASE(mForwardAlphaTestOnlyPS);
	SAFE_RELEASE(mForwardAlphaTestPS);
	SAFE_RELEASE(mForwardPS);
	SAFE_RELEASE(mClusteredPS);
	SAFE_RELEASE(mGBufferAlphaTestPS);
	SAFE_RELEASE(mGBufferPS);
	SAFE_RELEASE(mGeometryVS);
}

void App::OnD3D11ResizedSwapChain(ID3D11Device* d3dDevice,
	const DXGI_SURFACE_DESC* backBufferDesc)
{
	mGBufferWidth = backBufferDesc->Width;
	mGBufferHeight = backBufferDesc->Height;

	// Create/recreate any textures related to screen size
	mGBuffer.resize(0);
	mGBufferRTV.resize(0);
	mGBufferSRV.resize(0);
	mLitBufferPS = 0;
	mLitBufferCS = 0;
	mDeferredLightingAccumBuffer = 0;
	mDepthBuffer = 0;
	SAFE_RELEASE(mDepthBufferReadOnlyDSV);

	DXGI_SAMPLE_DESC sampleDesc;
	sampleDesc.Count = mMSAASamples;
	sampleDesc.Quality = 0;

	// standard depth/stencil buffer
	mDepthBuffer = shared_ptr<Depth2D>(new Depth2D(
		d3dDevice, mGBufferWidth, mGBufferHeight,
		D3D11_BIND_DEPTH_STENCIL | D3D11_BIND_SHADER_RESOURCE,
		sampleDesc,
		mMSAASamples > 1 // Include stencil if using MSAA
	));

	// read-only depth stencil view
	{
		D3D11_DEPTH_STENCIL_VIEW_DESC desc;
		mDepthBuffer->GetDepthStencil()->GetDesc(&desc);
		desc.Flags = D3D11_DSV_READ_ONLY_DEPTH;

		d3dDevice->CreateDepthStencilView(mDepthBuffer->GetTexture(), &desc, &mDepthBufferReadOnlyDSV);
	}

	// NOTE: The next set of buffers are not all needed at the same time... a given technique really only needs one of them.
	// We allocate them all up front for quick swapping between techniques and to keep the code as simple as possible.

	// lit buffers
	mLitBufferPS = shared_ptr<Texture2D>(new Texture2D(
		d3dDevice, mGBufferWidth, mGBufferHeight, DXGI_FORMAT_R16G16B16A16_FLOAT,
		D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE,
		sampleDesc));

	mLitBufferCS = shared_ptr<StructuredBuffer<FramebufferFlatElement>>(new StructuredBuffer<FramebufferFlatElement>(
		d3dDevice, mGBufferWidth * mGBufferHeight * mMSAASamples,
		D3D11_BIND_UNORDERED_ACCESS | D3D11_BIND_SHADER_RESOURCE));

	// deferred lighting accumulation buffer
	mDeferredLightingAccumBuffer = shared_ptr<Texture2D>(new Texture2D(
		d3dDevice, mGBufferWidth, mGBufferHeight, DXGI_FORMAT_R16G16B16A16_FLOAT,
		D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE,
		sampleDesc));

	// G-Buffer

	// normal_specular
	mGBuffer.push_back(shared_ptr<Texture2D>(new Texture2D(
		d3dDevice, mGBufferWidth, mGBufferHeight, DXGI_FORMAT_R16G16B16A16_FLOAT,
		D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE,
		sampleDesc)));

	// albedo
	mGBuffer.push_back(shared_ptr<Texture2D>(new Texture2D(
		d3dDevice, mGBufferWidth, mGBufferHeight, DXGI_FORMAT_R8G8B8A8_UNORM,
		D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE,
		sampleDesc)));

	// positionZgrad
	mGBuffer.push_back(shared_ptr<Texture2D>(new Texture2D(
		d3dDevice, mGBufferWidth, mGBufferHeight, DXGI_FORMAT_R16G16_FLOAT,
		D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE,
		sampleDesc)));

	// Set up GBuffer resource list
	mGBufferRTV.resize(mGBuffer.size(), 0);
	mGBufferSRV.resize(mGBuffer.size() + 1, 0);
	for (std::size_t i = 0; i < mGBuffer.size(); ++i)
	{
		mGBufferRTV[i] = mGBuffer[i]->GetRenderTarget();
		mGBufferSRV[i] = mGBuffer[i]->GetShaderResource();
	}
	// Depth buffer is the last SRV that we use for reading
	mGBufferSRV.back() = mDepthBuffer->GetShaderResource();

	// invalidate cache as it depends on framebuffer resolution for some reason
	mCachedCameraView._11 = 0;
}

void App::InitializeLightParameters(ID3D11Device* d3dDevice)
{
	mPointLightParameters.resize(MAX_LIGHTS);
	mLightInitialTransform.resize(MAX_LIGHTS);
	mPointLightPositionWorld.resize(MAX_LIGHTS);

	// Use a constant seed for consistency
	std::tr1::mt19937 rng(1337);

	std::tr1::uniform_real<float> radiusNormDist(0.0f, 1.0f);
	const float maxRadius = 100.0f;
	std::tr1::uniform_real<float> angleDist(0.0f, 2.0f * DirectX::XM_PI);
	std::tr1::uniform_real<float> heightDist(0.0f, 20.0f);
	std::tr1::uniform_real<float> animationSpeedDist(2.0f, 20.0f);
	std::tr1::uniform_int<int> animationDirection(0, 1);
	std::tr1::uniform_real<float> hueDist(0.0f, 1.0f);
	std::tr1::uniform_real<float> intensityDist(0.1f, 0.5f);
	std::tr1::uniform_real<float> attenuationDist(2.0f, 15.0f);
	const float attenuationStartFactor = 0.8f;

	for (unsigned int i = 0; i < MAX_LIGHTS; ++i)
	{
		PointLight& params = mPointLightParameters[i];
		PointLightInitTransform& init = mLightInitialTransform[i];

		init.radius = std::sqrt(radiusNormDist(rng)) * maxRadius;
		init.angle = angleDist(rng);
		init.height = heightDist(rng);
		// Normalize by arc length
		init.animationSpeed = (animationDirection(rng) * 2 - 1) * animationSpeedDist(rng) / init.radius;

		// HSL->RGB, vary light hue
		params.color = HueToRGB(hueDist(rng));
		params.color.x *= intensityDist(rng);
		params.color.y *= intensityDist(rng);
		params.color.z *= intensityDist(rng);
		params.attenuationEnd = attenuationDist(rng);
		params.attenuationBegin = attenuationStartFactor * params.attenuationEnd;
	}
}

void App::SetActiveLights(ID3D11Device* d3dDevice, double log2ActiveLights)
{
	int activeLights = int(pow(2, log2ActiveLights));
	mActiveLights = activeLights;

	delete mLightBuffer;
	mLightBuffer = new StructuredBuffer<PointLight>(d3dDevice, activeLights, D3D11_BIND_SHADER_RESOURCE, true);

	// Make sure all the active lights are set up
	Move(0.0f);
}

void App::Move(float elapsedTime)
{
	mTotalTime += elapsedTime;

	// Update positions of active lights
	for (unsigned int i = 0; i < mActiveLights; ++i)
	{
		const PointLightInitTransform& initTransform = mLightInitialTransform[i];
		float angle = initTransform.angle + mTotalTime * initTransform.animationSpeed;
		mPointLightPositionWorld[i] = DirectX::XMFLOAT3(
			initTransform.radius * std::cos(angle),
			initTransform.height,
			initTransform.radius * std::sin(angle));
	}

	// invalidate cache as it depends on light positions
	mCachedCameraView._11 = 0;
}

void App::Render(ID3D11DeviceContext* d3dDeviceContext,
	ID3D11RenderTargetView* backBuffer,
	CDXUTSDKMesh& mesh_opaque,
	CDXUTSDKMesh& mesh_alpha,
	ID3D11ShaderResourceView* skybox,
	const DirectX::XMFLOAT4X4A& worldMatrix,
	const CFirstPersonCamera* viewerCamera,
	const D3D11_VIEWPORT* viewport,
	const UIConstants* ui)
{
	DirectX::XMFLOAT4X4A cameraProj;
	DirectX::XMStoreFloat4x4A(&cameraProj, viewerCamera->GetProjMatrix());
	DirectX::XMFLOAT4X4A cameraView;
	DirectX::XMStoreFloat4x4A(&cameraView, viewerCamera->GetViewMatrix());

	DirectX::XMFLOAT4X4A cameraViewInv;
	DirectX::XMStoreFloat4x4A(&cameraViewInv, DirectX::XMMatrixInverse(NULL, viewerCamera->GetViewMatrix()));

	// Compute composite matrices
	DirectX::XMFLOAT4X4A cameraViewProj;
	DirectX::XMStoreFloat4x4(&cameraViewProj, DirectX::XMMatrixMultiply(DirectX::XMLoadFloat4x4A(&cameraView), DirectX::XMLoadFloat4x4A(&cameraProj)));

	DirectX::XMFLOAT4X4A cameraWorldViewProj;
	DirectX::XMStoreFloat4x4(&cameraWorldViewProj, DirectX::XMMatrixMultiply(DirectX::XMLoadFloat4x4A(&worldMatrix), DirectX::XMLoadFloat4x4A(&cameraViewProj)));

	// Fill in frame constants
	{
		D3D11_MAPPED_SUBRESOURCE mappedResource;
		d3dDeviceContext->Map(mPerFrameConstants, 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedResource);
		PerFrameConstants* constants = static_cast<PerFrameConstants*>(mappedResource.pData);

		constants->mCameraWorldViewProj = cameraWorldViewProj;
		DirectX::XMStoreFloat4x4(&constants->mCameraWorldView, DirectX::XMMatrixMultiply(DirectX::XMLoadFloat4x4A(&worldMatrix), DirectX::XMLoadFloat4x4A(&cameraView)));
		constants->mCameraViewProj = cameraViewProj;
		constants->mCameraProj = cameraProj;
		// NOTE: Complementary Z => swap near/far back
		constants->mCameraNearFar = DirectX::XMFLOAT4(viewerCamera->GetFarClip(), viewerCamera->GetNearClip(), 0.0f, 0.0f);

		constants->mFramebufferDimensionsX = mGBufferWidth;
		constants->mFramebufferDimensionsY = mGBufferHeight;
		constants->mFramebufferDimensionsZ = 0; // Unused
		constants->mFramebufferDimensionsW = 0; // Unused

		constants->mUI = *ui;

		d3dDeviceContext->Unmap(mPerFrameConstants, 0);
	}

	// Geometry phase
	if (mesh_opaque.IsLoaded())
	{
		//    mesh_opaque.ComputeInFrustumFlags(cameraWorldViewProj);
	}
	if (mesh_alpha.IsLoaded())
	{
		//    mesh_alpha.ComputeInFrustumFlags(cameraWorldViewProj);
	}

	// Setup lights
	ID3D11ShaderResourceView* lightBufferSRV = SetupLights(d3dDeviceContext, cameraView);

	// Forward rendering takes a different path here
	if (ui->lightCullTechnique == CULL_FORWARD_NONE)
	{
		RenderForward(d3dDeviceContext, mesh_opaque, mesh_alpha, lightBufferSRV, viewerCamera, viewport, ui, false);
	}
	else if (ui->lightCullTechnique == CULL_FORWARD_PREZ_NONE)
	{
		RenderForward(d3dDeviceContext, mesh_opaque, mesh_alpha, lightBufferSRV, viewerCamera, viewport, ui, true);
	}
	else if (ui->lightCullTechnique == CULL_CLUSTERED)
	{
		ID3D11ShaderResourceView* LightGridBufferSRV = SetupLightGrid(d3dDeviceContext, viewerCamera, ui);
		RenderClustered(d3dDeviceContext, mesh_opaque, mesh_alpha, lightBufferSRV, LightGridBufferSRV, viewerCamera, viewport, ui);
	}
	else
	{
		RenderGBuffer(d3dDeviceContext, mesh_opaque, mesh_alpha, viewerCamera, viewport, ui);
		ComputeLighting(d3dDeviceContext, lightBufferSRV, viewport, ui);
	}

	// Render skybox and tonemap
	RenderSkyboxAndToneMap(d3dDeviceContext, backBuffer, skybox,
		mDepthBuffer->GetShaderResource(), viewport, ui);
}

ID3D11ShaderResourceView* App::SetupLights(ID3D11DeviceContext* d3dDeviceContext,
	const DirectX::XMFLOAT4X4A& cameraView)
{
	// Transform light world positions into view space and store in our parameters array
	DirectX::XMVector3TransformCoordStream(
		&mPointLightParameters[0].positionView, sizeof(PointLight),
		&mPointLightPositionWorld[0], sizeof(DirectX::XMFLOAT3), mActiveLights, DirectX::XMLoadFloat4x4A(&cameraView));

	// Copy light list into shader buffer
	{
		PointLight* light = mLightBuffer->MapDiscard(d3dDeviceContext);
		for (unsigned int i = 0; i < mActiveLights; ++i)
		{
			light[i] = mPointLightParameters[i];
		}
		mLightBuffer->Unmap(d3dDeviceContext);
	}

	return mLightBuffer->GetShaderResource();
}

ID3D11ShaderResourceView* App::RenderForward(ID3D11DeviceContext* d3dDeviceContext,
	CDXUTSDKMesh& mesh_opaque,
	CDXUTSDKMesh& mesh_alpha,
	ID3D11ShaderResourceView* lightBufferSRV,
	const CFirstPersonCamera* viewerCamera,
	const D3D11_VIEWPORT* viewport,
	const UIConstants* ui,
	bool doPreZ)
{
	// Clear lit and depth buffer
	const float zeros[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
	d3dDeviceContext->ClearRenderTargetView(mLitBufferPS->GetRenderTarget(), zeros);
	// NOTE: Complementary Z buffer: clear to 0 (far)!
	d3dDeviceContext->ClearDepthStencilView(mDepthBuffer->GetDepthStencil(), D3D11_CLEAR_DEPTH, 0.0f, 0);

	d3dDeviceContext->IASetInputLayout(mMeshVertexLayout);

	d3dDeviceContext->VSSetConstantBuffers(0, 1, &mPerFrameConstants);
	d3dDeviceContext->VSSetShader(mGeometryVS, 0, 0);

	d3dDeviceContext->GSSetShader(0, 0, 0);

	d3dDeviceContext->RSSetViewports(1, viewport);

	d3dDeviceContext->PSSetConstantBuffers(0, 1, &mPerFrameConstants);
	d3dDeviceContext->PSSetShaderResources(5, 1, &lightBufferSRV);
	d3dDeviceContext->PSSetSamplers(0, 1, &mDiffuseSampler);
	// Diffuse texture set per-material by DXUT mesh routines

	d3dDeviceContext->OMSetDepthStencilState(mDepthState, 0);

	// Pre-Z pass if requested
	if (doPreZ)
	{
		d3dDeviceContext->OMSetRenderTargets(0, 0, mDepthBuffer->GetDepthStencil());

		// Render opaque geometry
		if (mesh_opaque.IsLoaded())
		{
			d3dDeviceContext->RSSetState(mRasterizerState);
			d3dDeviceContext->PSSetShader(0, 0, 0);
			mesh_opaque.Render(d3dDeviceContext, 0);
		}

		// Render alpha tested geometry
		if (mesh_alpha.IsLoaded())
		{
			d3dDeviceContext->RSSetState(mDoubleSidedRasterizerState);
			// NOTE: Use simplified alpha test shader that only clips
			d3dDeviceContext->PSSetShader(mForwardAlphaTestOnlyPS, 0, 0);
			mesh_alpha.Render(d3dDeviceContext, 0);
		}
	}

	// Set up render targets
	ID3D11RenderTargetView* renderTargets[1] = { mLitBufferPS->GetRenderTarget() };
	d3dDeviceContext->OMSetRenderTargets(1, renderTargets, mDepthBuffer->GetDepthStencil());
	d3dDeviceContext->OMSetBlendState(mGeometryBlendState, 0, 0xFFFFFFFF);

	// Render opaque geometry
	if (mesh_opaque.IsLoaded())
	{
		d3dDeviceContext->RSSetState(mRasterizerState);
		d3dDeviceContext->PSSetShader(mForwardPS, 0, 0);
		mesh_opaque.Render(d3dDeviceContext, 0);
	}

	// Render alpha tested geometry
	if (mesh_alpha.IsLoaded())
	{
		d3dDeviceContext->RSSetState(mDoubleSidedRasterizerState);
		d3dDeviceContext->PSSetShader(mForwardAlphaTestPS, 0, 0);
		mesh_alpha.Render(d3dDeviceContext, 0);
	}

	// Cleanup (aka make the runtime happy)
	d3dDeviceContext->OMSetRenderTargets(0, 0, 0);

	return mLitBufferPS->GetShaderResource();
}

UIConstants cachedUI;

ID3D11ShaderResourceView* App::SetupLightGrid(ID3D11DeviceContext* d3dDeviceContext, const CFirstPersonCamera* viewerCamera, const UIConstants* ui)
{
	int n = ui->clusteredGridScale;
	mLightGridBuilder.reset(LightGridDimensions(2 * n, n, 8 * n));

	// Uncomment this to enable grid caching
	// if (*viewerCamera->GetViewMatrix() != mCachedCameraView || memcmp(&cachedUI, ui, sizeof(UIConstants)) != 0)
	{
		DirectX::XMStoreFloat4x4(&mCachedCameraView, viewerCamera->GetViewMatrix());
		cachedUI = *ui;

		mLightGridBuilder.clearAllFragments();
		RasterizeLights(&mLightGridBuilder, viewerCamera, &mPointLightParameters[0], mActiveLights);

		LightGridEntry* gpuBuffer = mLightGridBuffer->MapDiscard(d3dDeviceContext);
		mLightGridBuilder.buildAndUpload(gpuBuffer, mLightGridBufferSize * 16);
		mLightGridBuffer->Unmap(d3dDeviceContext);
	}

	return mLightGridBuffer->GetShaderResource();
}

ID3D11ShaderResourceView* App::RenderClustered(ID3D11DeviceContext* d3dDeviceContext,
	CDXUTSDKMesh& mesh_opaque,
	CDXUTSDKMesh& mesh_alpha,
	ID3D11ShaderResourceView* lightBufferSRV,
	ID3D11ShaderResourceView* LightGridBufferSRV,
	const CFirstPersonCamera* viewerCamera,
	const D3D11_VIEWPORT* viewport,
	const UIConstants* ui)
{
	// Clear lit and depth buffer
	const float zeros[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
	d3dDeviceContext->ClearRenderTargetView(mLitBufferPS->GetRenderTarget(), zeros);
	// NOTE: Complementary Z buffer: clear to 0 (far)!
	d3dDeviceContext->ClearDepthStencilView(mDepthBuffer->GetDepthStencil(), D3D11_CLEAR_DEPTH, 0.0f, 0);

	d3dDeviceContext->IASetInputLayout(mMeshVertexLayout);

	d3dDeviceContext->VSSetConstantBuffers(0, 1, &mPerFrameConstants);
	d3dDeviceContext->VSSetShader(mGeometryVS, 0, 0);

	d3dDeviceContext->GSSetShader(0, 0, 0);

	d3dDeviceContext->RSSetViewports(1, viewport);

	d3dDeviceContext->PSSetConstantBuffers(0, 1, &mPerFrameConstants);
	d3dDeviceContext->PSSetShaderResources(5, 1, &lightBufferSRV);
	d3dDeviceContext->PSSetShaderResources(6, 1, &LightGridBufferSRV);
	d3dDeviceContext->PSSetSamplers(0, 1, &mDiffuseSampler);
	// Diffuse texture set per-material by DXUT mesh routines

	d3dDeviceContext->OMSetDepthStencilState(mDepthState, 0);

	// Set up render targets
	ID3D11RenderTargetView* renderTargets[1] = { mLitBufferPS->GetRenderTarget() };
	d3dDeviceContext->OMSetRenderTargets(1, renderTargets, mDepthBuffer->GetDepthStencil());
	d3dDeviceContext->OMSetBlendState(mGeometryBlendState, 0, 0xFFFFFFFF);

	// Render opaque geometry
	if (mesh_opaque.IsLoaded())
	{
		d3dDeviceContext->RSSetState(mRasterizerState);
		d3dDeviceContext->PSSetShader(mClusteredPS, 0, 0);
		mesh_opaque.Render(d3dDeviceContext, 0);
	}

	// Render alpha tested geometry
	if (mesh_alpha.IsLoaded())
	{
		assert(false && "not implemented");
	}

	// Cleanup (aka make the runtime happy)
	d3dDeviceContext->OMSetRenderTargets(0, 0, 0);

	return mLitBufferPS->GetShaderResource();
}

void App::RenderGBuffer(ID3D11DeviceContext* d3dDeviceContext,
	CDXUTSDKMesh& mesh_opaque,
	CDXUTSDKMesh& mesh_alpha,
	const CFirstPersonCamera* viewerCamera,
	const D3D11_VIEWPORT* viewport,
	const UIConstants* ui)
{
	// Clear GBuffer
	// NOTE: We actually only need to clear the depth buffer here since we replace unwritten (i.e. far plane) samples
	// with the skybox. We use the depth buffer to reconstruct position and only in-frustum positions are shaded.
	// NOTE: Complementary Z buffer: clear to 0 (far)!
	d3dDeviceContext->ClearDepthStencilView(mDepthBuffer->GetDepthStencil(), D3D11_CLEAR_DEPTH | D3D11_CLEAR_STENCIL, 0.0f, 0);

	d3dDeviceContext->IASetInputLayout(mMeshVertexLayout);

	d3dDeviceContext->VSSetConstantBuffers(0, 1, &mPerFrameConstants);
	d3dDeviceContext->VSSetShader(mGeometryVS, 0, 0);

	d3dDeviceContext->GSSetShader(0, 0, 0);

	d3dDeviceContext->RSSetViewports(1, viewport);

	d3dDeviceContext->PSSetConstantBuffers(0, 1, &mPerFrameConstants);
	d3dDeviceContext->PSSetSamplers(0, 1, &mDiffuseSampler);
	// Diffuse texture set per-material by DXUT mesh routines

	// Set up render GBuffer render targets
	d3dDeviceContext->OMSetDepthStencilState(mDepthState, 0);
	d3dDeviceContext->OMSetRenderTargets(static_cast<UINT>(mGBufferRTV.size()), &mGBufferRTV.front(), mDepthBuffer->GetDepthStencil());
	d3dDeviceContext->OMSetBlendState(mGeometryBlendState, 0, 0xFFFFFFFF);

	// Render opaque geometry
	if (mesh_opaque.IsLoaded())
	{
		d3dDeviceContext->RSSetState(mRasterizerState);
		d3dDeviceContext->PSSetShader(mGBufferPS, 0, 0);
		mesh_opaque.Render(d3dDeviceContext, 0);
	}

	// Render alpha tested geometry
	if (mesh_alpha.IsLoaded())
	{
		d3dDeviceContext->RSSetState(mDoubleSidedRasterizerState);
		d3dDeviceContext->PSSetShader(mGBufferAlphaTestPS, 0, 0);
		mesh_alpha.Render(d3dDeviceContext, 0);
	}

	// Cleanup (aka make the runtime happy)
	d3dDeviceContext->OMSetRenderTargets(0, 0, 0);
}

void App::CullComputeShaderTile(ID3D11DeviceContext* d3dDeviceContext,
	ID3D11ShaderResourceView* lightBufferSRV)
{
	// No need to clear, we write all pixels

	// Compute shader setup (always does all the lights at once)
	d3dDeviceContext->CSSetConstantBuffers(0, 1, &mPerFrameConstants);
	d3dDeviceContext->CSSetShaderResources(0, static_cast<UINT>(mGBufferSRV.size()), &mGBufferSRV.front());
	d3dDeviceContext->CSSetShaderResources(5, 1, &lightBufferSRV);

	ID3D11UnorderedAccessView* litBufferUAV = mLitBufferCS->GetUnorderedAccess();
	d3dDeviceContext->CSSetUnorderedAccessViews(0, 1, &litBufferUAV, 0);
	d3dDeviceContext->CSSetShader(mComputeShaderTileCS, 0, 0);

	// Dispatch
	unsigned int dispatchWidth = (mGBufferWidth + COMPUTE_SHADER_TILE_GROUP_DIM - 1) / COMPUTE_SHADER_TILE_GROUP_DIM;
	unsigned int dispatchHeight = (mGBufferHeight + COMPUTE_SHADER_TILE_GROUP_DIM - 1) / COMPUTE_SHADER_TILE_GROUP_DIM;
	d3dDeviceContext->Dispatch(dispatchWidth, dispatchHeight, 1);
}

void App::CullQuad(ID3D11DeviceContext* d3dDeviceContext,
	ID3D11ShaderResourceView* lightBufferSRV,
	const D3D11_VIEWPORT* viewport,
	const UIConstants* ui)
{
	bool deferredLighting = (ui->lightCullTechnique == CULL_QUAD_DEFERRED_LIGHTING);
	std::tr1::shared_ptr<Texture2D>& accumulateBuffer = deferredLighting ? mDeferredLightingAccumBuffer : mLitBufferPS;

	// Clear
	const float zeros[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
	d3dDeviceContext->ClearRenderTargetView(accumulateBuffer->GetRenderTarget(), zeros);

	if (mMSAASamples > 1)
	{
		// Full screen triangle setup
		d3dDeviceContext->IASetInputLayout(0);
		d3dDeviceContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
		d3dDeviceContext->IASetVertexBuffers(0, 0, 0, 0, 0);

		d3dDeviceContext->VSSetShader(mFullScreenTriangleVS, 0, 0);
		d3dDeviceContext->GSSetShader(0, 0, 0);

		d3dDeviceContext->RSSetState(mRasterizerState);
		d3dDeviceContext->RSSetViewports(1, viewport);

		d3dDeviceContext->PSSetConstantBuffers(0, 1, &mPerFrameConstants);
		d3dDeviceContext->PSSetShaderResources(0, static_cast<UINT>(mGBufferSRV.size()), &mGBufferSRV.front());
		d3dDeviceContext->PSSetShaderResources(5, 1, &lightBufferSRV);

		// Set stencil mask for samples that require per-sample shading
		d3dDeviceContext->PSSetShader(mRequiresPerSampleShadingPS, 0, 0);
		d3dDeviceContext->OMSetDepthStencilState(mWriteStencilState, 1);
		d3dDeviceContext->OMSetRenderTargets(0, 0, mDepthBufferReadOnlyDSV);
		d3dDeviceContext->Draw(3, 0);
	}

	// Point primitives expanded into quads in the geometry shader
	d3dDeviceContext->IASetInputLayout(0);
	d3dDeviceContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
	d3dDeviceContext->IASetVertexBuffers(0, 0, 0, 0, 0);

	d3dDeviceContext->VSSetConstantBuffers(0, 1, &mPerFrameConstants);
	d3dDeviceContext->VSSetShaderResources(5, 1, &lightBufferSRV);
	d3dDeviceContext->VSSetShader(mGPUQuadVS, 0, 0);

	d3dDeviceContext->GSSetShader(mGPUQuadGS, 0, 0);

	d3dDeviceContext->RSSetState(mRasterizerState);
	d3dDeviceContext->RSSetViewports(1, viewport);

	d3dDeviceContext->PSSetConstantBuffers(0, 1, &mPerFrameConstants);
	d3dDeviceContext->PSSetShaderResources(0, static_cast<UINT>(mGBufferSRV.size()), &mGBufferSRV.front());
	d3dDeviceContext->PSSetShaderResources(5, 1, &lightBufferSRV);

	// Additively blend into lit buffer
	ID3D11RenderTargetView* renderTargets[1] = { accumulateBuffer->GetRenderTarget() };
	// Use depth buffer for culling but no writes (use the read-only DSV)
	d3dDeviceContext->OMSetRenderTargets(1, renderTargets, mDepthBufferReadOnlyDSV);
	d3dDeviceContext->OMSetBlendState(mLightingBlendState, 0, 0xFFFFFFFF);

	// Dispatch one point per light

	// Do pixel frequency shading
	d3dDeviceContext->PSSetShader(deferredLighting ? mGPUQuadDLPS : mGPUQuadPS, 0, 0);
	d3dDeviceContext->OMSetDepthStencilState(mEqualStencilState, 0);
	d3dDeviceContext->Draw(mActiveLights, 0);

	if (mMSAASamples > 1)
	{
		// Do sample frequency shading
		d3dDeviceContext->PSSetShader(deferredLighting ? mGPUQuadDLPerSamplePS : mGPUQuadPerSamplePS, 0, 0);
		d3dDeviceContext->OMSetDepthStencilState(mEqualStencilState, 1);
		d3dDeviceContext->Draw(mActiveLights, 0);
	}

	if (deferredLighting)
	{
		// Final screen-space pass to combine diffuse and specular
		d3dDeviceContext->IASetInputLayout(0);
		d3dDeviceContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
		d3dDeviceContext->IASetVertexBuffers(0, 0, 0, 0, 0);

		d3dDeviceContext->VSSetShader(mFullScreenTriangleVS, 0, 0);
		d3dDeviceContext->GSSetShader(0, 0, 0);

		ID3D11RenderTargetView* resolveRenderTargets[1] = { mLitBufferPS->GetRenderTarget() };
		d3dDeviceContext->OMSetRenderTargets(1, resolveRenderTargets, mDepthBufferReadOnlyDSV);
		d3dDeviceContext->OMSetBlendState(mGeometryBlendState, 0, 0xFFFFFFFF);

		ID3D11ShaderResourceView* accumulateBufferSRV = accumulateBuffer->GetShaderResource();
		d3dDeviceContext->PSSetShaderResources(7, 1, &accumulateBufferSRV);

		// Do pixel frequency resolve
		d3dDeviceContext->PSSetShader(mGPUQuadDLResolvePS, 0, 0);
		d3dDeviceContext->OMSetDepthStencilState(mEqualStencilState, 0);
		d3dDeviceContext->Draw(3, 0);

		if (mMSAASamples > 1)
		{
			// Do sample frequency resolve
			d3dDeviceContext->PSSetShader(mGPUQuadDLResolvePerSamplePS, 0, 0);
			d3dDeviceContext->OMSetDepthStencilState(mEqualStencilState, 1);
			d3dDeviceContext->Draw(3, 0);
		}
	}
}

void App::CullDeferredNone(ID3D11DeviceContext* d3dDeviceContext,
	ID3D11ShaderResourceView* lightBufferSRV,
	const D3D11_VIEWPORT* viewport)
{
	// Clear
	const float zeros[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
	d3dDeviceContext->ClearRenderTargetView(mLitBufferPS->GetRenderTarget(), zeros);

	// Full screen triangle setup
	d3dDeviceContext->IASetInputLayout(0);
	d3dDeviceContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
	d3dDeviceContext->IASetVertexBuffers(0, 0, 0, 0, 0);

	d3dDeviceContext->VSSetShader(mFullScreenTriangleVS, 0, 0);
	d3dDeviceContext->GSSetShader(0, 0, 0);

	d3dDeviceContext->RSSetState(mRasterizerState);
	d3dDeviceContext->RSSetViewports(1, viewport);

	d3dDeviceContext->PSSetConstantBuffers(0, 1, &mPerFrameConstants);
	d3dDeviceContext->PSSetShaderResources(0, static_cast<UINT>(mGBufferSRV.size()), &mGBufferSRV.front());
	d3dDeviceContext->PSSetShaderResources(5, 1, &lightBufferSRV);

	if (mMSAASamples > 1)
	{
		// Set stencil mask for samples that require per-sample shading
		d3dDeviceContext->PSSetShader(mRequiresPerSampleShadingPS, 0, 0);
		d3dDeviceContext->OMSetDepthStencilState(mWriteStencilState, 1);
		d3dDeviceContext->OMSetRenderTargets(0, 0, mDepthBufferReadOnlyDSV);
		d3dDeviceContext->Draw(3, 0);
	}

	// Additively blend into back buffer
	ID3D11RenderTargetView* renderTargets[1] = { mLitBufferPS->GetRenderTarget() };
	d3dDeviceContext->OMSetRenderTargets(1, renderTargets, mDepthBufferReadOnlyDSV);
	d3dDeviceContext->OMSetBlendState(mLightingBlendState, 0, 0xFFFFFFFF);

	// Do pixel frequency shading
	d3dDeviceContext->PSSetShader(mBasicLoopPS, 0, 0);
	d3dDeviceContext->OMSetDepthStencilState(mEqualStencilState, 0);
	d3dDeviceContext->Draw(3, 0);

	if (mMSAASamples > 1)
	{
		// Do sample frequency shading
		d3dDeviceContext->PSSetShader(mBasicLoopPerSamplePS, 0, 0);
		d3dDeviceContext->OMSetDepthStencilState(mEqualStencilState, 1);
		d3dDeviceContext->Draw(3, 0);
	}
}

void App::ComputeLighting(ID3D11DeviceContext* d3dDeviceContext,
	ID3D11ShaderResourceView* lightBufferSRV,
	const D3D11_VIEWPORT* viewport,
	const UIConstants* ui)
{
	switch (ui->lightCullTechnique)
	{
	case CULL_COMPUTE_SHADER_TILE:
	{
		CullComputeShaderTile(d3dDeviceContext, lightBufferSRV);
	}
	break;

	case CULL_QUAD:
	case CULL_QUAD_DEFERRED_LIGHTING:
	{
		CullQuad(d3dDeviceContext,
			lightBufferSRV,
			viewport,
			ui);
	}

	break;

	case CULL_DEFERRED_NONE:
	{
		CullDeferredNone(d3dDeviceContext,
			lightBufferSRV,
			viewport);
	}
	break;

	}; // switch

	// Cleanup (aka make the runtime happy)
	d3dDeviceContext->VSSetShader(0, 0, 0);
	d3dDeviceContext->GSSetShader(0, 0, 0);
	d3dDeviceContext->PSSetShader(0, 0, 0);
	d3dDeviceContext->OMSetRenderTargets(0, 0, 0);
	ID3D11ShaderResourceView* nullSRV[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	d3dDeviceContext->VSSetShaderResources(0, 8, nullSRV);
	d3dDeviceContext->PSSetShaderResources(0, 8, nullSRV);
	d3dDeviceContext->CSSetShaderResources(0, 8, nullSRV);
	ID3D11UnorderedAccessView* nullUAV[1] = { 0 };
	d3dDeviceContext->CSSetUnorderedAccessViews(0, 1, nullUAV, 0);
}

void App::RenderSkyboxAndToneMap(ID3D11DeviceContext* d3dDeviceContext,
	ID3D11RenderTargetView* backBuffer,
	ID3D11ShaderResourceView* skybox,
	ID3D11ShaderResourceView* depthSRV,
	const D3D11_VIEWPORT* viewport,
	const UIConstants* ui)
{
	D3D11_VIEWPORT skyboxViewport(*viewport);
	skyboxViewport.MinDepth = 1.0f;
	skyboxViewport.MaxDepth = 1.0f;

	d3dDeviceContext->IASetInputLayout(mMeshVertexLayout);

	d3dDeviceContext->VSSetConstantBuffers(0, 1, &mPerFrameConstants);
	d3dDeviceContext->VSSetShader(mSkyboxVS, 0, 0);

	d3dDeviceContext->RSSetState(mDoubleSidedRasterizerState);
	d3dDeviceContext->RSSetViewports(1, &skyboxViewport);

	d3dDeviceContext->PSSetConstantBuffers(0, 1, &mPerFrameConstants);
	d3dDeviceContext->PSSetSamplers(0, 1, &mDiffuseSampler);
	d3dDeviceContext->PSSetShader(mSkyboxPS, 0, 0);

	d3dDeviceContext->PSSetShaderResources(5, 1, &skybox);
	d3dDeviceContext->PSSetShaderResources(6, 1, &depthSRV);

	// Bind the appropriate lit buffer depending on the technique
	ID3D11ShaderResourceView* litViews[2] = { 0, 0 };
	switch (ui->lightCullTechnique)
	{
		// Compute-shader based techniques use the flattened MSAA buffer
	case CULL_COMPUTE_SHADER_TILE:
		litViews[1] = mLitBufferCS->GetShaderResource();
		break;
	default:
		litViews[0] = mLitBufferPS->GetShaderResource();
		break;
	}
	d3dDeviceContext->PSSetShaderResources(7, 2, litViews);

	d3dDeviceContext->OMSetRenderTargets(1, &backBuffer, 0);
	d3dDeviceContext->OMSetBlendState(mGeometryBlendState, 0, 0xFFFFFFFF);

	mSkyboxMesh.Render(d3dDeviceContext);

	// Cleanup (aka make the runtime happy)
	d3dDeviceContext->OMSetRenderTargets(0, 0, 0);
	ID3D11ShaderResourceView* nullViews[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	d3dDeviceContext->PSSetShaderResources(0, 10, nullViews);
}
