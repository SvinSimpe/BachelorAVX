cbuffer PER_FRAME_DATA
{
	float4x4	viewMatrix;
	float4x4	projMatrix;
	float3		eyePosition;
};

cbuffer PER_OBJECT_DATA
{
	float4x4 worldMatrix;
};


struct VertexInput
{
	float3 position	: POSITION;
};

struct VertexOutput
{
	float4 position	: SV_POSITION;
};

// VertexShader just passes through the input vertex without performing
// any modifications on it. Transformations occur in GeometryShader.
VertexInput VS_main( VertexInput input )
{
	return input;
}

[maxvertexcount(4)]
void GS_main( point VertexInput input[1], inout TriangleStream<VertexOutput> tStream )
{
	VertexOutput temp = (VertexOutput)0;

	float3 camToParticle = normalize( input[0].position.xyz - eyePosition );
	float3 right		 = normalize( cross( float3( 0.0f, 1.0f, 0.0f ), camToParticle ) );
	float3 up			 = normalize( cross( camToParticle, right ) );

	float size = 0.5f;
	
	for( int i = 0; i < 4; i++ )
	{
		temp.position.xyz = input[0].position;
		if( i == 0 )
			temp.position += float4( -right - up, 0.0f ) * size;

		else if( i == 1 )
			temp.position += float4( -right + up, 0.0f ) * size;

		else if( i == 2 )
			temp.position += float4(  right - up, 0.0f ) * size;
		
		else
			temp.position += float4(  right + up, 0.0f ) * size;

		temp.position = mul( temp.position, worldMatrix );
		temp.position = mul( temp.position, viewMatrix );
		temp.position = mul( temp.position, projMatrix );

		tStream.Append( temp );
	}
}

float4 PS_main( VertexOutput input ) : SV_Target
{
	// Bright White Snow, HO HO HO
	return float4( 1.0f, 1.0f, 1.0f, 1.0f );
}