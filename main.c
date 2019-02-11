#define _CRT_SECURE_NO_WARNINGS

#include "lcms2.h"
#include "md5.h"

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// SMPTE ST.2084: https://ieeexplore.ieee.org/servlet/opac?punumber=7291450

static const float PQ_C1 = 0.8359375;       // 3424.0 / 4096.0
static const float PQ_C2 = 18.8515625;      // 2413.0 / 4096.0 * 32.0
static const float PQ_C3 = 18.6875;         // 2392.0 / 4096.0 * 32.0
static const float PQ_M1 = 0.1593017578125; // 2610.0 / 4096.0 / 4.0
static const float PQ_M2 = 78.84375;        // 2523.0 / 4096.0 * 128.0

// SMPTE ST.2084: Equation 4.1
// L = ( (max(N^(1/m2) - c1, 0)) / (c2 - c3*N^(1/m2)) )^(1/m1)
static float PQ_EOTF(float N)
{
    float N1m2 = powf(N, 1 / PQ_M2);
    float N1m2c1 = N1m2 - PQ_C1;
    if (N1m2c1 < 0.0f)
        N1m2c1 = 0.0f;
    float c2c3N1m2 = PQ_C2 - (PQ_C3 * N1m2);
    return powf(N1m2c1 / c2c3N1m2, 1 / PQ_M1);
}

void makePQLUT(float * out)
{
    const float c1 = 0.8359375f;       // 3424 / 4096
    const float c2 = 18.8515625f;      // 2413 / 4096 * 32
    const float c3 = 18.6875f;         // 2392 / 4096 * 32
    const float m1 = 0.1593017578125f; // 2610 / 4096 / 4
    const float m2 = 78.84375f;        // 2523 / 4096 * 128

    for (int i = 0; i < 4096; ++i) {
        float src = (float)i / 4095.0f;
        out[i] = PQ_EOTF(src);
    }
}

int main(int argc, char * argv[])
{
    float primaries[8];
    primaries[0] = 0.64f;
    primaries[1] = 0.33f;
    primaries[2] = 0.30f;
    primaries[3] = 0.60f;
    primaries[4] = 0.15f;
    primaries[5] = 0.06f;
    primaries[6] = 0.3127f;
    primaries[7] = 0.3290f;

    float gamma = 2.2f;
    float maxLuminance = 10000;

    cmsContext lcms = cmsCreateContext(NULL, NULL);

    cmsToneCurve * curves[3];
    cmsCIExyYTRIPLE dstPrimaries;
    cmsCIExyY dstWhitePoint;
    cmsCIEXYZ lumi;

    dstPrimaries.Red.x = primaries[0];
    dstPrimaries.Red.y = primaries[1];
    dstPrimaries.Red.Y = 0.0f; // unused
    dstPrimaries.Green.x = primaries[2];
    dstPrimaries.Green.y = primaries[3];
    dstPrimaries.Green.Y = 0.0f; // unused
    dstPrimaries.Blue.x = primaries[4];
    dstPrimaries.Blue.y = primaries[5];
    dstPrimaries.Blue.Y = 0.0f; // unused
    dstWhitePoint.x = primaries[6];
    dstWhitePoint.y = primaries[7];
    dstWhitePoint.Y = 1.0f;

    float * toneCurve = (float *)malloc(sizeof(float) * 4096);
    makePQLUT(toneCurve);

    curves[0] = cmsBuildTabulatedToneCurveFloat(lcms, 4096, toneCurve); //cmsBuildGamma(lcms, gamma);
    curves[1] = curves[0];
    curves[2] = curves[0];

    cmsHPROFILE profile = cmsCreateRGBProfileTHR(lcms, &dstWhitePoint, &dstPrimaries, curves);
    cmsFreeToneCurve(curves[0]);

    lumi.X = 0.0f;
    lumi.Y = maxLuminance;
    lumi.Z = 0.0f;
    cmsWriteTag(profile, cmsSigLuminanceTag, &lumi);

    cmsUInt32Number bytesNeeded;
    cmsSaveProfileToMem(profile, NULL, &bytesNeeded);
    char * raw = malloc(bytesNeeded);
    cmsSaveProfileToMem(profile, raw, &bytesNeeded);

    FILE * out = fopen("out.icc", "wb");
    fwrite(raw, 1, bytesNeeded, out);
    fclose(out);

    cmsHPROFILE rereadProfile = cmsOpenProfileFromMemTHR(lcms, raw, bytesNeeded);
    int rawCurveSize = cmsReadRawTag(rereadProfile, cmsSigRedTRCTag, NULL, 0);
    if (rawCurveSize > 0) {
        uint8_t * rawCurve = malloc(rawCurveSize);
        cmsReadRawTag(rereadProfile, cmsSigRedTRCTag, rawCurve, rawCurveSize);

        FILE * outCurve = fopen("pqCurve.bin", "wb");
        fwrite(rawCurve, 1, rawCurveSize, outCurve);
        fclose(outCurve);

        cmsToneCurve * rereadCurve = (cmsToneCurve *)cmsReadTag(rereadProfile, cmsSigRedTRCTag);
        float pqGamma = (float)cmsEstimateGamma(rereadCurve, 1.0f);
        printf("Estimated PQ gamma: %f\n", pqGamma);

        {
            uint8_t signature[16];

            MD5_CTX ctx;
            MD5_Init(&ctx);
            MD5_Update(&ctx, rawCurve, (unsigned long)rawCurveSize);
            MD5_Final(signature, &ctx);

            uint8_t * s = signature;
            printf("MD5: %x%x%x%x%x%x%x%x%x%x%x%x%x%x%x%x\n",
                s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15]);
        }
    }

    return 0;
}
