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

#ifndef _APP_H_
#define _APP_H_ 1

#define WIN32_LEAN_AND_MEAN 1
#define NOMINMAX 1
#include <DXUT.h>
#include <DXUTcamera.h>
#include <SDKMesh.h>
#include <DirectXMath.h>
#include "Texture2D.h"
#include "Buffer.h"
#include "LightGrid.h"
#include <vector>
#include <memory>

enum LightCullTechnique {
	CULL_FORWARD_NONE = 0,
	CULL_FORWARD_PREZ_NONE,
	CULL_DEFERRED_NONE,
	CULL_QUAD,
	CULL_QUAD_DEFERRED_LIGHTING,
	CULL_COMPUTE_SHADER_TILE,
	CULL_CLUSTERED,
};

// NOTE: Must match shader equivalent structure
__declspec(align(16))
struct UIConstants
{
	unsigned int lightingOnly;
	unsigned int faceNormals;
	unsigned int visualizeLightCount;
	unsigned int visualizePerSampleShading;
	unsigned int lightCullTechnique;

	unsigned int clusteredGridScale;
};

// NOTE: Must match shader equivalent structure
struct PointLight
{
	DirectX::XMFLOAT3 positionView;
	float attenuationBegin;
	DirectX::XMFLOAT3 color;
	float attenuationEnd;
};


// Host-side world-space initial transform data
struct PointLightInitTransform
{
	// Cylindrical coordinates
	float radius;
	float angle;
	float height;
	float animationSpeed;
};

// Flat framebuffer RGBA16-encoded
struct FramebufferFlatElement
{
	unsigned int rb;
	unsigned int ga;
};


class App
{
public:
	App(ID3D11Device* d3dDevice, double log2ActiveLights, unsigned int msaaSamples);

	~App();

	void OnD3D11ResizedSwapChain(ID3D11Device* d3dDevice,
		const DXGI_SURFACE_DESC* backBufferDesc);

	void Move(float elapsedTime);

	void Render(ID3D11DeviceContext* d3dDeviceContext,
		ID3D11RenderTargetView* backBuffer,
		CDXUTSDKMesh& mesh_opaque,
		CDXUTSDKMesh& mesh_alpha,
		ID3D11ShaderResourceView* skybox,
		const DirectX::XMFLOAT4X4A& worldMatrix,
		const CFirstPersonCamera* viewerCamera,
		const D3D11_VIEWPORT* viewport,
		const UIConstants* ui);

	void SetActiveLights(ID3D11Device* d3dDevice, double log2ActiveLights);
	unsigned int GetActiveLights() const { return mActiveLights; }

private:
	void InitializeLightParameters(ID3D11Device* d3dDevice);

	// Notes: 
	// - Most of these functions should all be called after initializing per frame/pass constants, etc.
	//   as the shaders that they invoke bind those constant buffers.

	// Set up shader light buffer
	ID3D11ShaderResourceView* SetupLights(ID3D11DeviceContext* d3dDeviceContext,
		const DirectX::XMFLOAT4X4A& cameraView);

	std::vector<uint32_t> CullLightFrustrum(const CFirstPersonCamera* viewerCamera, int x, int y, int z, int grid[3]);

	ID3D11ShaderResourceView* SetupLightGrid(ID3D11DeviceContext* d3dDeviceContext,
		const CFirstPersonCamera* viewerCamera,
		const UIConstants* ui);

	// Forward rendering of geometry into
	ID3D11ShaderResourceView* RenderForward(ID3D11DeviceContext* d3dDeviceContext,
		CDXUTSDKMesh& mesh_opaque,
		CDXUTSDKMesh& mesh_alpha,
		ID3D11ShaderResourceView* lightBufferSRV,
		const CFirstPersonCamera* viewerCamera,
		const D3D11_VIEWPORT* viewport,
		const UIConstants* ui,
		bool doPreZ);

	// Forward rendering of geometry into
	ID3D11ShaderResourceView* RenderClustered(ID3D11DeviceContext* d3dDeviceContext,
		CDXUTSDKMesh& mesh_opaque,
		CDXUTSDKMesh& mesh_alpha,
		ID3D11ShaderResourceView* lightBufferSRV,
		ID3D11ShaderResourceView* LightGridBufferSRV,
		const CFirstPersonCamera* viewerCamera,
		const D3D11_VIEWPORT* viewport,
		const UIConstants* ui);

	// Draws geometry into G-buffer
	void RenderGBuffer(ID3D11DeviceContext* d3dDeviceContext,
		CDXUTSDKMesh& mesh_opaque,
		CDXUTSDKMesh& mesh_alpha,
		const CFirstPersonCamera* viewerCamera,
		const D3D11_VIEWPORT* viewport,
		const UIConstants* ui);

	// Handles skybox, tone mapping, etc
	void RenderSkyboxAndToneMap(ID3D11DeviceContext* d3dDeviceContext,
		ID3D11RenderTargetView* backBuffer,
		ID3D11ShaderResourceView* skybox,
		ID3D11ShaderResourceView* depthSRV,
		const D3D11_VIEWPORT* viewport,
		const UIConstants* ui);

	void App::CullComputeShaderTile(ID3D11DeviceContext* d3dDeviceContext,
		ID3D11ShaderResourceView* lightBufferSRV);

	void App::CullQuad(ID3D11DeviceContext* d3dDeviceContext,
		ID3D11ShaderResourceView* lightBufferSRV,
		const D3D11_VIEWPORT* viewport,
		const UIConstants* ui);

	void App::CullDeferredNone(ID3D11DeviceContext* d3dDeviceContext,
		ID3D11ShaderResourceView* lightBufferSRV,
		const D3D11_VIEWPORT* viewport);


	void ComputeLighting(ID3D11DeviceContext* d3dDeviceContext,
		ID3D11ShaderResourceView* lightBufferSRV,
		const D3D11_VIEWPORT* viewport,
		const UIConstants* ui);

	unsigned int mMSAASamples;
	float mTotalTime;

	ID3D11InputLayout* mMeshVertexLayout;

	ID3D11VertexShader* mGeometryVS;

	ID3D11PixelShader* mGBufferPS;
	ID3D11PixelShader* mGBufferAlphaTestPS;

	ID3D11PixelShader* mForwardPS;
	ID3D11PixelShader* mForwardAlphaTestPS;
	ID3D11PixelShader* mForwardAlphaTestOnlyPS;

	ID3D11PixelShader* mClusteredPS;

	CDXUTSDKMesh mSkyboxMesh;
	ID3D11VertexShader* mSkyboxVS;
	ID3D11PixelShader* mSkyboxPS;

	ID3D11VertexShader* mFullScreenTriangleVS;

	ID3D11PixelShader* mRequiresPerSampleShadingPS;

	ID3D11PixelShader* mBasicLoopPS;
	ID3D11PixelShader* mBasicLoopPerSamplePS;

	ID3D11ComputeShader* mComputeShaderTileCS;

	ID3D11VertexShader* mGPUQuadVS;
	ID3D11GeometryShader* mGPUQuadGS;
	ID3D11PixelShader* mGPUQuadPS;
	ID3D11PixelShader* mGPUQuadPerSamplePS;

	ID3D11PixelShader* mGPUQuadDLPS;
	ID3D11PixelShader* mGPUQuadDLPerSamplePS;

	ID3D11PixelShader* mGPUQuadDLResolvePS;
	ID3D11PixelShader* mGPUQuadDLResolvePerSamplePS;

	ID3D11Buffer* mPerFrameConstants;

	ID3D11RasterizerState* mRasterizerState;
	ID3D11RasterizerState* mDoubleSidedRasterizerState;

	ID3D11DepthStencilState* mDepthState;
	ID3D11DepthStencilState* mWriteStencilState;
	ID3D11DepthStencilState* mEqualStencilState;

	ID3D11BlendState* mGeometryBlendState;
	ID3D11BlendState* mLightingBlendState;

	ID3D11SamplerState* mDiffuseSampler;

	std::vector< std::shared_ptr<Texture2D> > mGBuffer;
	// Handy cache of list of RT pointers for G-buffer
	std::vector<ID3D11RenderTargetView*> mGBufferRTV;
	// Handy cache of list of SRV pointers for the G-buffer
	std::vector<ID3D11ShaderResourceView*> mGBufferSRV;
	unsigned int mGBufferWidth;
	unsigned int mGBufferHeight;

	// We use a different lit buffer (different bind flags and MSAA handling) depending on whether we
	// write to it from the pixel shader (render target) or compute shader (UAV)
	std::shared_ptr<Texture2D> mLitBufferPS;
	std::shared_ptr<StructuredBuffer<FramebufferFlatElement> > mLitBufferCS;

	// A temporary accumulation buffer used for deferred lighting
	std::shared_ptr<Texture2D> mDeferredLightingAccumBuffer;

	std::shared_ptr<Depth2D> mDepthBuffer;
	// We also need a read-only depth stencil view for techniques that read the G-buffer while also using Z-culling
	ID3D11DepthStencilView* mDepthBufferReadOnlyDSV;

	// Lighting state
	unsigned int mActiveLights;
	std::vector<PointLightInitTransform> mLightInitialTransform;
	std::vector<PointLight> mPointLightParameters;
	std::vector<DirectX::XMFLOAT3> mPointLightPositionWorld;

	StructuredBuffer<PointLight>* mLightBuffer;

	// Clustered Shading state
	LightGridBuilder mLightGridBuilder;
	int mLightGridBufferSize;
	StructuredBuffer<LightGridEntry>* mLightGridBuffer;

	DirectX::XMFLOAT4X4 mCachedCameraView;
};

#endif