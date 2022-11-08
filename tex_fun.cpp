/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

GzColor	*image=NULL;
int xs, ys;
int reset = 1;

inline int ARRAY(int x, int y) { return (x + y * xs); }
/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd);
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
  // u,v => [0, 1]
  if (u > 1.0f) {
      u = 1.0f;
  }
  if (u < 0) {
      u = 0;
  }

  if (v > 1.0f) {
      v = 1.0f;
  }
  if (v < 0) {
      v = 0;
  }
/* determine texture cell corner values and perform bilinear interpolation */
  /* scaling u, v to framebuffer requirement */
  float factorU = (float)(xs - 1);
  float factorV = (float)(ys - 1);
  float scaledU = u * factorU;
  float scaledV = v * factorV;

  /* Corner values of u, v */
  GzColor A, B, C, D;

  int lowerU = floor(scaledU);
  int upperU = ceil(scaledU);
  int lowerV = floor(scaledV);
  int upperV = ceil(scaledV);


  int pointA = ARRAY(lowerU, lowerV);
  int pointB = ARRAY(upperU, lowerV);
  int pointC = ARRAY(upperU, upperV);
  int pointD = ARRAY(lowerU, upperV);

  for (int i = 0; i < 3; i++) {
      A[i] = image[pointA][i];
  }

  for (int i = 0; i < 3; i++) {
      B[i] = image[pointB][i];
  }

  for (int i = 0; i < 3; i++) {
      C[i] = image[pointC][i];
  }

  for (int i = 0; i < 3; i++) {
      D[i] = image[pointD][i];
  }

  /* Bilinear interpolation */
  /* using the equation:
  * color = (s * t) * C + (1 - s) * t * D + s * (1 - t) * B + (1 - s) * (1 - t) * A
  */
  float s = scaledU - (float)floor(scaledU);
  float t = scaledV - (float)floor(scaledV);
  /* set color to interpolated GzColor value and return */
  for (int i = 0; i < 3; i++) {
      color[i] = s * t * C[i] + (1.0f - s) * t * D[i] + s * (1.0f - t) * B[i] + (1.0f - s) * (1.0f - t) * A[i];
  }
	return GZ_SUCCESS;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
    if (u > 1.0f) {
        u = 1.0f;
    }
    if (u < 0) {
        u = 0;
    }
    if (v > 1.0f) {
        v = 1.0f;
    }
    if (v < 0) {
        v = 0;
    }
    int BLACK = 0;
    int WHITE = 4095;
    int interval = 6;
    int ceiledU = ceil(u * interval);
    int ceiledV = ceil(v * interval);


    if (ceiledU % 2 == 0 && ceiledV % 2 == 0) {
        color[RED] = WHITE;
        color[GREEN] = WHITE;
        color[BLUE] = WHITE;
    }
    else if (ceiledU % 2 != 0 && ceiledV % 2 != 0) {
        color[RED] = WHITE;
        color[GREEN] = WHITE;
        color[BLUE] = WHITE;
    }
    else {
        color[RED] = WHITE;
        color[GREEN] = BLACK;
        color[BLUE] = BLACK;
    }


    /* sine function procedural */
    //color[RED] = 2 * 50 * u * v;
    //color[GREEN] = 2 * sin(50 * u * v);
    //color[BLUE] = 2 * cos(50 * u * v);
     

	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

