// XNA math code ported from http://blogs.msdn.com/b/chuckw/archive/2012/07/28/spherical-harmonics-math.aspx
// Original code under MS-PL

//-------------------------------------------------------------------------------------
// DirectXSH.cpp -- C++ Spherical Harmonics Math Library
//
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//  
// Copyright (c) Microsoft Corporation. All rights reserved.
//-------------------------------------------------------------------------------------

#define XM_PI               3.141592654
#define XM_PIDIV2 (XM_PI * 0.5)
const int XM_SH_MINORDER = 2;
const int XM_SH_MAXORDER = 6;

#define REAL float
#define const
#define CONSTANT(x) (float(x))

//const float fExtraNormFac[XM_SH_MAXORDER] = { 2.0*sqrt(XM_PI), 2.0/3.0*sqrt(3.0*XM_PI), 2.0f/5.0*sqrt(5.0*XM_PI), 2.0/7.0*sqrt(7.0*XM_PI), 2.0/3.0*sqrt(XM_PI), 2.0/11.0*sqrt(11.0f*XM_PI) };
float fExtraNormFac(int v)
{
	float fv = float(v * 2 + 1);
	return 2.0*fv*sqrt(fv * XM_PI);
}

// routine generated programmatically for evaluating SH basis for degree 1
// inputs (x,y,z) are a point on the sphere (i.e., must be unit length)
// output is vector b with SH basis evaluated at (x,y,z).
//
void sh_eval_basis_1(REAL x,REAL y,REAL z,out REAL b[4])
{
	/* m=0 */
	
	// l=0
	const REAL p_0_0 = CONSTANT(0.282094791773878140);
	b[  0] = p_0_0; // l=0,m=0
	// l=1
	const REAL p_1_0 = CONSTANT(0.488602511902919920)*z;
	b[  2] = p_1_0; // l=1,m=0
	
	
	/* m=1 */
	
	const REAL s1 = y;
	const REAL c1 = x;
	
	// l=1
	const REAL p_1_1 = CONSTANT(-0.488602511902919920);
	b[  1] = p_1_1*s1; // l=1,m=-1
	b[  3] = p_1_1*c1; // l=1,m=+1
}

// routine generated programmatically for evaluating SH basis for degree 2
// inputs (x,y,z) are a point on the sphere (i.e., must be unit length)
// output is vector b with SH basis evaluated at (x,y,z).
//
void sh_eval_basis_2(REAL x,REAL y,REAL z,out REAL b[9])
{
	const REAL z2 = z*z;
	
	
	/* m=0 */
	
	// l=0
	const REAL p_0_0 = CONSTANT(0.282094791773878140);
	b[  0] = p_0_0; // l=0,m=0
	// l=1
	const REAL p_1_0 = CONSTANT(0.488602511902919920)*z;
	b[  2] = p_1_0; // l=1,m=0
	// l=2
	const REAL p_2_0 = CONSTANT(0.946174695757560080)*z2 + CONSTANT(-0.315391565252520050);
	b[  6] = p_2_0; // l=2,m=0
	
	
	/* m=1 */
	
	const REAL s1 = y;
	const REAL c1 = x;
	
	// l=1
	const REAL p_1_1 = CONSTANT(-0.488602511902919920);
	b[  1] = p_1_1*s1; // l=1,m=-1
	b[  3] = p_1_1*c1; // l=1,m=+1
	// l=2
	const REAL p_2_1 = CONSTANT(-1.092548430592079200)*z;
	b[  5] = p_2_1*s1; // l=2,m=-1
	b[  7] = p_2_1*c1; // l=2,m=+1
	
	
	/* m=2 */
	
	const REAL s2 = x*s1 + y*c1;
	const REAL c2 = x*c1 - y*s1;
	
	// l=2
	const REAL p_2_2 = CONSTANT(0.546274215296039590);
	b[  4] = p_2_2*s2; // l=2,m=-2
	b[  8] = p_2_2*c2; // l=2,m=+2
}

// input pF only consists of Yl0 values, normalizes coefficients for directional
// lights.
float CosWtInt(const int order)
{
	const float fCW0 = 0.25;
	const float fCW1 = 0.5;
	const float fCW2 = 5.0/16.0;
	//const float fCW3 = 0.0;
	const float fCW4 = -3.0/32.0;
	//const float fCW5 = 0.0;
	
	// order has to be at least linear...
	
	float fRet = fCW0 + fCW1;
	
	if (order > 2) fRet += fCW2;
	if (order > 4) fRet += fCW4;
	
	// odd degrees >= 3 evaluate to zero integrated against cosine...
	
	return fRet;
}

// computes the integral of a constant function over a solid angular
// extent.  No error checking - only used internaly.  This function
// only returns the Yl0 coefficients, since the rest are zero for
// circularly symmetric functions.
const float ComputeCapInt_t1 = sqrt(0.3141593E1);
const float ComputeCapInt_t5 = sqrt(3.0);
const float ComputeCapInt_t11 = sqrt(5.0);
const float ComputeCapInt_t18 = sqrt(7.0);
const float ComputeCapInt_t32 = sqrt(11.0);

void ComputeCapInt(const int order, float angle, out float pR[6])
{
	const float t2 = cos(angle);
	const float t3 = ComputeCapInt_t1*t2;
	const float t7 = sin(angle);
	const float t8 = t7*t7;
	
	for(int i = 0 ; i < 6 ; i++)
		pR[i] = 0.0;
	
	pR[0] = -t3+ComputeCapInt_t1;
	pR[1] = ComputeCapInt_t5*ComputeCapInt_t1*t8/2.0;
	
	if (order > 2)
	{
		const float t13 = t2*t2;
		
		pR[2] = -ComputeCapInt_t11*ComputeCapInt_t1*t2*(t13-1.0)/2.0;
		if (order > 3)
		{
			const float t19 = ComputeCapInt_t18*ComputeCapInt_t1;
			const float t20 = t13*t13;
			
			pR[3] = -5.0/8.0*t19*t20+3.0/4.0*t19*t13-t19/8.0;
			if (order > 4)
			{
				pR[4] = -3.0/8.0*t3*(7.0*t20-10.0*t13+3.0);
				if (order > 5)
				{
					const float t33 = ComputeCapInt_t32*ComputeCapInt_t1;
					pR[5] = -21.0/16.0*t33*t20*t13+35.0/16.0*t33*t20-15.0/16.0*t33*t13+t33/16.0;
				}
			}
		}
	}
}

void XMSHEvalDirectionalLight(vec3 dir, out float result[2 * 2])
{
    sh_eval_basis_1(dir.x, dir.y, dir.z, result); // evaluate the BF in this direction...

    // now compute "normalization" and scale vector for each valid spectral band
    const float fNorm = XM_PI / CosWtInt(2);
    
    for( int i=0; i < 2 * 2; ++i)
        result[i] *= fNorm;
}

void XMSHEvalDirectionalLight(vec3 dir, out float result[3 * 3])
{
    sh_eval_basis_2(dir.x, dir.y, dir.z, result); // evaluate the BF in this direction...

    // now compute "normalization" and scale vector for each valid spectral band
    const float fNorm = XM_PI / CosWtInt(3);
    
    for( int i=0; i < 3 * 3; ++i)
        result[i] *= fNorm;
}

void XMSHEvalSphericalLight(vec3 pos, float radius, out float result[2 * 2])
{
    const float fDist = length( pos );

    // WARNING: fDist should not be < radius - otherwise light contains origin

    //const float fSinConeAngle = (fDist <= radius) ? 0.99999f : radius/fDist;
    const float fConeAngle = (fDist <= radius) ? (XM_PIDIV2) : asin(radius/fDist);

    vec3 dir = normalize( pos );

    float fTmpL0[ XM_SH_MAXORDER ];

    const float fNewNorm = 1.0;///(fSinConeAngle*fSinConeAngle); 

    ComputeCapInt(2,fConeAngle,fTmpL0);

    const float fX = dir.x;
    const float fY = dir.y;
    const float fZ = dir.z;

	sh_eval_basis_1(fX,fY,fZ,result);
	
    const float fValUse0 = fTmpL0[0]*fNewNorm*fExtraNormFac(0);
	result[0] *= fValUse0;
    const float fValUse1 = fTmpL0[1]*fNewNorm*fExtraNormFac(1);
	for(int i = 1 ; i < 4 ; i++)
		result[i] *= fValUse1;
}

void XMSHEvalSphericalLight(vec3 pos, float radius, out float result[3 * 3])
{
    const float fDist = length( pos );

    // WARNING: fDist should not be < radius - otherwise light contains origin

    //const float fSinConeAngle = (fDist <= radius) ? 0.99999f : radius/fDist;
    const float fConeAngle = (fDist <= radius) ? (XM_PIDIV2) : asin(radius/fDist);

    vec3 dir = normalize( pos );

    float fTmpL0[ XM_SH_MAXORDER ];

    const float fNewNorm = 1.0;///(fSinConeAngle*fSinConeAngle); 

    ComputeCapInt(3,fConeAngle,fTmpL0);

    const float fX = dir.x;
    const float fY = dir.y;
    const float fZ = dir.z;

	sh_eval_basis_2(fX,fY,fZ,result);
	
    const float fValUse0 = fTmpL0[0]*fNewNorm*fExtraNormFac(0);
	result[0] *= fValUse0;
    const float fValUse1 = fTmpL0[1]*fNewNorm*fExtraNormFac(1);
	for(int i = 1 ; i < 4 ; i++)
		result[i] *= fValUse1;
    const float fValUse2 = fTmpL0[2]*fNewNorm*fExtraNormFac(2);
	for(int i = 4 ; i < 9 ; i++)
		result[i] *= fValUse2;
}

///////////////////////////////////////////////////////////////////////
// Example code for applying SH
///////////////////////////////////////////////////////////////////////

float Eval2(float sh_l[2 * 2], float sh_n[2 * 2])
{
	// Lambertian diffuse per-band convolution coeffs
	float sh_c[2];
	sh_c[0] = 1.0;
	sh_c[1] = 2.0/3.0;
	
	// Evaluate exit radiance in direction of normal
	float c = 0.0;
	/*
	for (int l = 0; l < 2; l++)
		for (int m = -l; m <= l; m++, i++)
			c += sh_c[l]*sh_n[i]*sh_l[i];
	*/
	c = sh_c[0] * (sh_n[0] * sh_l[0]);
	for(int i = 1 ; i < 4 ; i++)
		c += sh_c[1] * sh_n[i] * sh_l[i];
	return c;
}

float Eval3(float sh_l[3 * 3], float sh_n[3 * 3])
{
	// Lambertian diffuse per-band convolution coeffs
	float sh_c[3];
	sh_c[0] = 1.0;
	sh_c[1] = 2.0/3.0;
	sh_c[2] = 1.0/4.0;
	
	// Evaluate exit radiance in direction of normal
	float c = 0.0;
	c = sh_c[0] * (sh_n[0] * sh_l[0]);
	for(int i = 1 ; i < 4 ; i++)
		c += sh_c[1] * sh_n[i] * sh_l[i];
	for(int i = 4 ; i < 9 ; i++)
		c += sh_c[2] * sh_n[i] * sh_l[i];
	return c;
}

///////////////////////////////////////////////////////////////////////
// Example code to calculate and preview SH
///////////////////////////////////////////////////////////////////////
vec4 GetSpherePosInBox(vec2 pos, vec2 center, float size)
{
	vec2 uv = (pos - center) / size;
	float offset = uv.x * uv.x + uv.y * uv.y;
	if(offset > 1.0)
		return vec4(0);

	float z = sqrt(1.0 - offset);
	return vec4(uv, z, 1.0);
}

const vec3 up = vec3(0.0, 1.0, 0.0);
float DirectionalOrder2(vec3 spherePos)
{
	float lightSH[2 * 2], normalSH[2 * 2];
	XMSHEvalDirectionalLight(up, lightSH);
	sh_eval_basis_1(spherePos.x, spherePos.y, spherePos.z, normalSH);
	
	return clamp(Eval2(lightSH, normalSH), 0.0, 1.0);
}

float DirectionalOrder3(vec3 spherePos)
{
	float lightSH[3 * 3], normalSH[3 * 3];
	XMSHEvalDirectionalLight(up, lightSH);
	sh_eval_basis_2(spherePos.x, spherePos.y, spherePos.z, normalSH);
	
	return clamp(Eval3(lightSH, normalSH), 0.0, 1.0);
}

float PointOrder2(vec3 spherePos)
{
	float lightSH[2 * 2], normalSH[2 * 2];
	XMSHEvalSphericalLight(4.0 * up, 1.0, lightSH);
	sh_eval_basis_1(spherePos.x, spherePos.y, spherePos.z, normalSH);
	
	return clamp(Eval2(lightSH, normalSH), 0.0, 1.0);
}

float PointOrder3(vec3 spherePos)
{
	float lightSH[3 * 3], normalSH[3 * 3];
	XMSHEvalSphericalLight(4.0 * up, 1.0, lightSH);
	sh_eval_basis_2(spherePos.x, spherePos.y, spherePos.z, normalSH);
	
	return clamp(Eval3(lightSH, normalSH), 0.0, 1.0);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 center = iResolution.xy * 0.5;
	float minAxis = 0.5 * min(iResolution.x, iResolution.y);
	bool right = fragCoord.x > center.x;
	bool top = fragCoord.y > center.y;
	vec2 gridCenter = iResolution.xy * 0.25;
	
	if(top)
		gridCenter.y += iResolution.y * 0.5;
	if(right)
		gridCenter.x += iResolution.x * 0.5;
	
	vec4 spherePos = GetSpherePosInBox(fragCoord.xy, gridCenter, minAxis * 0.5);
	float result = 0.0;
	if(top)
	{
		if(right)
			result = DirectionalOrder3(spherePos.xyz);
		else
			result = DirectionalOrder2(spherePos.xyz);
		
		if(fract(iTime) > 0.5)
		{
			result -= clamp(dot(spherePos.xyz, up), 0.0, 1.0);
			result = abs(result);// * 10.0;
		}	
	}
	else
	{
		if(right)
			result = PointOrder3(spherePos.xyz);
		else
			result = PointOrder2(spherePos.xyz);
	}
	vec4 bgColor = 0.5 * vec4(float(top), float(right), float(!top&&!right), 0.0);
	fragColor = mix(bgColor, vec4(result), spherePos.w);
}