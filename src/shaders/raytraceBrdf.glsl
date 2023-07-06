
#include "shared_structs.h"

#define PI 3.141592

vec3 EvalBrdf(vec3 N, vec3 L, vec3 V, Material mat)
{
    vec3 Kd = mat.diffuse;
    vec3 Ks = mat.specular;
    const float alpha = mat.shininess;

    vec3 H = normalize(L + V);
    float LH = dot(L, H);

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

    // G factor with V
    float mV = dot(H, V);
    float NV = dot(N, V);
    float tan_square_theta_v = (1.0 -  NV * NV) /  NV * NV;
    float GV;
    if(NV > 1.0 || sqrt(tan_square_theta_v) == 0)
    {
        GV = 1.0;
    }
    else
    {
        int x;
        if (mV / NV > 0) x = 1;
        else x = 0;
        GV = x * 2 / (1.0 + sqrt(1 + alpha_square * tan_square_theta_v));
    }

    // G factor with L
    float mL = LH; // dot(L, H);
    float NL = dot(N, L);
    float tan_square_theta_l = (1.0 -  NL * NL) /  NL * NL;
    float GL;
    if(NL > 1.0 || sqrt(tan_square_theta_l) == 0)
    {
        GL = 1.0;
    }
    else
    {
        int x;
        if (mL / NL > 0) x = 1;
        else x = 0;
        GL = x * 2 / (1.0 + sqrt(1 + alpha_square * tan_square_theta_l));
    }
    float G = GV * GL;

    return max(NL, 0.0) * ((Kd / PI) + ( (D * G * F) / (4 * abs(NL) * abs(NV)) ) );
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