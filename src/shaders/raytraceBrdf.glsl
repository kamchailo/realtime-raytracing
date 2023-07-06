
#include "shared_structs.h"

#define PI 3.141592

vec3 EvalBrdf(vec3 N, vec3 L, vec3 V, Material mat)
{
    vec3 Kd = mat.diffuse;
    vec3 Ks = mat.specular;
    const float alpha = mat.shininess;

    vec3 H = normalize(L + V);
    float LH = max(dot(L, H), 0.0);

    // L = Wi
    // V = Wo
    // m = H

    // F factor
    vec3 F = Ks + (vec3(1.0) - Ks) * pow((1 - LH), 5);

    // D factor
    float mN = dot(H, N);
    float tan_square_theta_m = (1.0 - mN * mN) / mN * mN;
    float alpha_square = alpha * alpha;
    float D = clamp(mN, 0.0, 1.0) * alpha_square
            / (PI * pow(mN, 4) * pow(alpha_square + tan_square_theta_m, 2));

    // G factor
    float mV = dot(H, V);
    float NV = dot(N, V);
    float tan_square_theta_v = (1.0 -  NV * NV) /  NV * NV;
    float G;
    if(NV > 1.0 || sqrt(tan_square_theta_v)==0)
    {
        G = 1.0;
    }
    else
    {
        G = clamp(mV/NV, 0.0, 1.0) * 2 / (1.0 + sqrt(1 + alpha_square * tan_square_theta_v));
    }

    // vec3 Ambient = pcRay.tempAmbient.xyz;
    // vec3 Ii = pcRay.tempLightInt.xyz;
    // return only BRDF
    float NL = dot(N, L);
    return max(NL, 0.0) * ((Kd / PI) + ( (D * G * F) / (4 * abs(NL) * abs(NV)) ) );
    
}

vec3 EvalBrdf2(vec3 N, vec3 L, vec3 V, Material mat)
{
    vec3 Kd = mat.diffuse;
    vec3 Ks = mat.specular;
    const float alpha = mat.shininess;
    vec3 H = normalize(L+V); // same as (Wo+Wi)||Wo+Wi||
    // L = wi
    // V = wo
    // H = m

    // D Term
    float dotmN = dot(H, N);
    float tan_theta_m_squared = (1.0 - dotmN * dotmN) / (dotmN * dotmN);
    float num =pow(alpha,2);
    float denom = PI * pow(dotmN,4) * pow((pow(alpha,2) + tan_theta_m_squared),2);
    float D = clamp(dotmN, 0.0, 1.0) * num/denom;

    // F Term 
    float dotLH = max(dot(L,H),0.0);
    vec3 F = Ks + (vec3(1.0) - Ks) * pow((1-dotLH),5);

    // G Term
    float dotvN = dot(V, N);
    float dotvm = dot(V, H);
    float tan_theta_v_squared = (1.0 - dotvN * dotvN) / (dotvN * dotvN);
    float G=0;
    if(dotvN > 1.0 || sqrt(tan_theta_v_squared)==0)
    {
        G = 1.0;
    }
    else
    {
        G = clamp(dotvm/dotvN, 0.0, 1.0) * 2 / (1 + sqrt(1+pow(alpha,2)* tan_theta_v_squared));
    }
    return ( max( dot(N,L), 0.0) * ( Kd/PI + (D*F*G) / ( 4 * abs(dot(V,N)) * abs(dot(L,N)) ) ) );//
}

vec3 SampleLobe(vec3 A, float c, float phi)
{
    float s = sqrt(1 - c * c); // sin angle

    // Create vector K centered around Z-axis and rotate to A-axis
    vec3 K = vec3( s * cos(phi), s * sin(phi), c);
    if( abs(A.z - 1.0) < 1e-3)
    {
        return K; // A == Z so no rotation
    }
    if( abs(A.z + 1.0) < 1e-3)
    {
        return vec3(K.x, -K.y, -K.z); // A == -Z so rotate 180 around X axis
    }

    // A = normalize(A); // Not needed if A is unit length
    vec3 B = normalize( vec3(-A.y, A.x, 0.0)); // B = Z x A
    vec3 C = cross(A, B);

    /*
        return vec3(1.0,0.0,0.0);
    */
	return K.x * B + K.y * C + K.z * A;
}

vec3 SampleBrdf(inout uint seed, in vec3 N)
{
	return SampleLobe(N, sqrt(rnd(seed)), 2 * PI * rnd(seed));
}

float PdfBrdf(vec3 N, vec3 Wi)
{
	return abs(dot(N, Wi)) / PI;
}