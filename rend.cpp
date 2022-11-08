/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include    "unordered_map"

#define PI (float) 3.14159265358979323846
#define RGB 3
#define LIMIT 4095


/* helper functions definitions */
void detailShade(GzCoord norm, GzCoord color, GzLight* lights, GzLight ambientLight, GzCoord camera, float spec, GzColor ks, GzColor kd, GzColor ka);
void textureShadePhong(GzCoord norm, GzColor color, GzLight* lights, GzCoord camera, GzColor ks, GzColor textureColor, GzLight ambientLight, float spec);
void textureShadeGouraud(GzCoord norm, GzColor color, GzLight* lights, GzCoord camera, GzColor textureColor, GzLight ambientLight, float spec);
float getMagnitude(float* vector);
float helperDotProduct(float* hor, float* vert);
GzIntensity boundCheck(GzIntensity pixelIntensity);
GzIntensity bitShiftingToChar(GzIntensity color);
float cosine(float degree);
float sine(float degree);
void dotProduct(GzMatrix matVert, GzMatrix transformationStack[], int level);
void getCrossProduct(float* v1, float* v2, float* result);
void makePureRotationMatrix(GzMatrix regMatrix, GzMatrix pureMatrix);
void normalize(GzCoord norm, GzCoord normalized);

/* vertex class for easier computation */
class Vertex {
public:
	GzCoord normal;
	GzColor color;
	GzCoord vertex;
	GzTextureIndex texture;

	void Vertex::setValues(GzCoord v) {
		vertex[X] = v[X];
		vertex[Y] = v[Y];
		vertex[Z] = v[Z];
	}
	void Vertex::setValuesColor(GzColor c, GzCoord v, GzTextureIndex uv) {
		color[RED] = c[RED];
		color[GREEN] = c[GREEN];
		color[BLUE] = c[BLUE];

		vertex[X] = v[X];
		vertex[Y] = v[Y];
		vertex[Z] = v[Z];

		texture[U] = uv[U];
		texture[V] = uv[V];
	}
	void Vertex::setValuesNormal(GzCoord norm, GzCoord v, GzTextureIndex uv) {
		normalize(norm, normal);

		vertex[X] = v[X];
		vertex[Y] = v[Y];
		vertex[Z] = v[Z];

		texture[U] = uv[U];
		texture[V] = uv[V];
	}

};

/* helper functions that need the custom class */
void sortVertices(Vertex vertices[3]);
void drawTriangleGouraudShading(GzRender* r, Vertex v1, Vertex v2, Vertex v3);
void drawTrianglePhongShading(GzRender* r, Vertex v1, Vertex v2, Vertex v3, GzLight* lights, GzLight ambient, GzCoord camera, float spec, GzColor ks, GzColor kd, GzColor ka);
void warp(Vertex* vert);
void unwarp(float z, GzTextureIndex warped, GzTextureIndex unwarpedUV);


/* Differential Digital Analyzer (for scan-line rasterization) - I created 3 different types depending on the shading style.
TODO: Create an interface for the DDA to support all the shading styles 
*/
class DDA {
public:
	/* geometry parameters */
	GzCoord start;
	GzCoord end;
	GzCoord current;

	/* normal parameters */
	GzCoord startNormal;
	GzCoord endNormal;
	GzCoord currentNormal;

	/* color parameters */
	GzColor startColor;
	GzColor endColor;
	GzColor currentColor;

	/*texture parameters*/
	GzTextureIndex startUV;
	GzTextureIndex endUV;
	GzTextureIndex currentUV;

	/* geometry interpolators */
	float slopeX;
	float slopeZ;

	/* normal interpolators */
	float normalXSlope;
	float normalYSlope;
	float normalZSlope;

	/* tesxture interpolators */
	float slopeU;
	float slopeV;

	/* color interpolators */
	float redSlope;
	float greenSlope;
	float blueSlope;

	/* regular DDA */
	void DDA::setValues(Vertex begin, Vertex finish) {
		start[0] = begin.vertex[0];
		start[1] = begin.vertex[1];
		start[2] = begin.vertex[2];

		end[0] = finish.vertex[0];
		end[1] = finish.vertex[1];
		end[2] = finish.vertex[2];


		slopeX = (end[0] - start[0]) / (end[1] - start[1]);
		slopeZ = (end[2] - start[2]) / (end[1] - start[1]);


		current[0] = start[0];
		current[1] = start[1];
		current[2] = start[2];

	}
	/* DDA for Phong shading (normal interpolation) */
	void DDA::setValuesNormal(Vertex begin, Vertex finish) {
		normalize(begin.normal, startNormal);
		normalize(finish.normal, endNormal);

		start[0] = begin.vertex[0];
		start[1] = begin.vertex[1];
		start[2] = begin.vertex[2];

		end[0] = finish.vertex[0];
		end[1] = finish.vertex[1];
		end[2] = finish.vertex[2];

		startUV[U] = begin.texture[U];
		startUV[V] = begin.texture[V];

		endUV[U] = finish.texture[U];
		endUV[V] = finish.texture[V];


		slopeU = (endUV[U] - startUV[U]) / (end[1] - start[1]);
		slopeV = (endUV[V] - startUV[V]) / (end[1] - start[1]);


		slopeX = (end[0] - start[0]) / (end[1] - start[1]);
		slopeZ = (end[2] - start[2]) / (end[1] - start[1]);
		normalXSlope = (endNormal[0] - startNormal[0]) / (end[1] - start[1]);
		normalYSlope = (endNormal[1] - startNormal[1]) / (end[1] - start[1]);
		normalZSlope = (endNormal[2] - startNormal[2]) / (end[1] - start[1]);


		current[0] = start[0];
		current[1] = start[1];
		current[2] = start[2];

		currentNormal[0] = startNormal[0];
		currentNormal[1] = startNormal[1];
		currentNormal[2] = startNormal[2];

		currentUV[U] = startUV[U];
		currentUV[V] = startUV[V];

	}
	/* DDA for Gouraud shading (color interpolation )*/
	void DDA::setValuesColor(Vertex begin, Vertex finish) {
		start[0] = begin.vertex[0];
		start[1] = begin.vertex[1];
		start[2] = begin.vertex[2];

		end[0] = finish.vertex[0];
		end[1] = finish.vertex[1];
		end[2] = finish.vertex[2];


		startColor[RED] = begin.color[RED];
		startColor[GREEN] = begin.color[GREEN];
		startColor[BLUE] = begin.color[BLUE];

		endColor[RED] = finish.color[RED];
		endColor[GREEN] = finish.color[GREEN];
		endColor[BLUE] = finish.color[BLUE];


		slopeX = (end[0] - start[0]) / (end[1] - start[1]);
		slopeZ = (end[2] - start[2]) / (end[1] - start[1]);
		redSlope = (endColor[RED] - startColor[RED]) / (end[1] - start[1]);
		greenSlope = (endColor[GREEN] - startColor[GREEN]) / (end[1] - start[1]);
		blueSlope = (endColor[BLUE] - startColor[BLUE]) / (end[1] - start[1]);

		startUV[U] = begin.texture[U];
		startUV[V] = begin.texture[V];

		endUV[U] = finish.texture[U];
		endUV[V] = finish.texture[V];


		slopeU = (endUV[U] - startUV[U]) / (end[1] - start[1]);
		slopeV = (endUV[V] - startUV[V]) / (end[1] - start[1]);

		currentUV[U] = startUV[U];
		currentUV[V] = startUV[V];


		current[X] = start[X];
		current[Y] = start[Y];
		current[Z] = start[Z];

		currentColor[RED] = startColor[RED];
		currentColor[GREEN] = startColor[GREEN];
		currentColor[BLUE] = startColor[BLUE];

	}

};


/* helper fucntions implementations */

/* to get the magnitude or length of a vector */
float getMagnitude(float* vector) {
	float sumSoFar = 0.0;
	for (int i = 0; i < 3; i++) {
		sumSoFar += pow(vector[i], 2);
	}
	float magnitude = sqrt(sumSoFar);
	return magnitude;
}

/* to get the similarity between 2 vectors */
float helperDotProduct(float* hor, float* vert) {
	float sumSoFar = 0.0;
	for (int i = 0; i < 3; i++) {
		sumSoFar += (hor[i] * vert[i]);
	}
	return sumSoFar;
}
/* helper function to shade textures */
void textureShadePhong(GzCoord norm, GzColor color, GzLight* lights, GzCoord camera, GzColor ks, GzColor textureColor, GzLight ambientLight, float spec) {
	GzColor specColor = { 0, 0, 0 };
	GzColor diffusedColor = { 0, 0, 0 };
	GzCoord R;
	

	for (int i = 0; i < 3; i++) {
		float nDotL = helperDotProduct(norm, lights[i].direction);
		float nDotE = helperDotProduct(norm, camera);
		/* flip normal */
		if (nDotE < 0 && nDotL < 0) {
			/* flip normal */
			GzCoord flippedNorm;
			for (int j = 0; j < 3; j++) {
				flippedNorm[j] = (-1.0f * norm[j]);
			}
			nDotL = helperDotProduct(flippedNorm, lights[i].direction);
			nDotE = helperDotProduct(flippedNorm, camera);

			/* calculate R */
			for (int j = 0; j < 3; j++) {
				R[j] = (flippedNorm[j] * 2 * nDotL) - lights[i].direction[j];
			}

			/* R dot E */
			float rDotE = helperDotProduct(R, camera);
			if (rDotE > 1) {
				rDotE = 1;
			}
			else if (rDotE < 0) {
				rDotE = 0;
			}
			float resultOfSpecPower = pow(rDotE, spec);

			for (int k = 0; k < 3; k++) {
				specColor[k] += (resultOfSpecPower * lights[i].color[k]);
				diffusedColor[k] += (nDotL * lights[i].color[k]);
			}

		}
		else if (nDotL > 0 && nDotE > 0) {

			/* calculate R */
			for (int j = 0; j < 3; j++) {
				R[j] = (norm[j] * 2 * nDotL) - lights[i].direction[j];
			}

			/* R dot E */
			float rDotE = helperDotProduct(R, camera);
			if (rDotE > 1) {
				rDotE = 1;
			}
			else if (rDotE < 0) {
				rDotE = 0;
			}
			float resultOfSpecPower = pow(rDotE, spec);
			for (int k = 0; k < 3; k++) {
				specColor[k] += (resultOfSpecPower * lights[i].color[k]);
				diffusedColor[k] += (nDotL * lights[i].color[k]);
			}
		}
	}
	/* calculate color */
	float c;
	for (int i = 0; i < 3; i++) {
		c = (ks[i] * specColor[i]) + (textureColor[i] * diffusedColor[i]) + (textureColor[i] * ambientLight.color[i]);
		if (c > 1) {
			color[i] = 1.0f;
		}
		else if (c < 0) {
			color[i] = 0;
		}
		else {
			color[i] = c;
		}
	}
}
/* helper function for gouraud shading with textures */
void textureShadeGouraud(GzCoord norm, GzColor color, GzLight* lights, GzCoord camera,  GzColor textureColor, GzLight ambientLight, float spec) {
	GzColor specColor = { 0, 0, 0 };
	GzColor diffusedColor = { 0, 0, 0 };
	GzCoord R;


	for (int i = 0; i < 3; i++) {
		float nDotL = helperDotProduct(norm, lights[i].direction);
		float nDotE = helperDotProduct(norm, camera);
		/* flip normal */
		if (nDotE < 0 && nDotL < 0) {
			/* flip normal */
			GzCoord flippedNorm;
			for (int j = 0; j < 3; j++) {
				flippedNorm[j] = (-1.0f * norm[j]);
			}
			nDotL = helperDotProduct(flippedNorm, lights[i].direction);
			nDotE = helperDotProduct(flippedNorm, camera);

			/* calculate R */
			for (int j = 0; j < 3; j++) {
				R[j] = (flippedNorm[j] * 2 * nDotL) - lights[i].direction[j];
			}

			/* R dot E */
			float rDotE = helperDotProduct(R, camera);
			if (rDotE > 1) {
				rDotE = 1;
			}
			else if (rDotE < 0) {
				rDotE = 0;
			}
			float resultOfSpecPower = pow(rDotE, spec);

			for (int k = 0; k < 3; k++) {
				specColor[k] += (resultOfSpecPower * lights[i].color[k]);
				diffusedColor[k] += (nDotL * lights[i].color[k]);
			}

		}
		else if (nDotL > 0 && nDotE > 0) {

			/* calculate R */
			for (int j = 0; j < 3; j++) {
				R[j] = (norm[j] * 2 * nDotL) - lights[i].direction[j];
			}

			/* R dot E */
			float rDotE = helperDotProduct(R, camera);
			if (rDotE > 1) {
				rDotE = 1;
			}
			else if (rDotE < 0) {
				rDotE = 0;
			}
			float resultOfSpecPower = pow(rDotE, spec);
			for (int k = 0; k < 3; k++) {
				specColor[k] += (resultOfSpecPower * lights[i].color[k]);
				diffusedColor[k] += (nDotL * lights[i].color[k]);
			}
		}
	}
	/* calculate color */
	float c;
	for (int i = 0; i < 3; i++) {
		c = (textureColor[i] * specColor[i]) + (textureColor[i] * diffusedColor[i]) + (textureColor[i] * ambientLight.color[i]);
		if (c > 1) {
			color[i] = 1.0f;
		}
		else if (c < 0) {
			color[i] = 0;
		}
		else {
			color[i] = c;
		}
	}
}

/* helper fucntion to shade (phong) */
void detailShade(GzCoord norm, GzCoord color, GzLight* lights, GzLight ambientLight, GzCoord camera, float spec, GzColor ks, GzColor kd, GzColor ka)
{
	GzColor specColor = { 0, 0, 0 };
	GzColor diffusedColor = { 0, 0, 0 };
	GzCoord R;


	for (int i = 0; i < 3; i++) {
		float nDotL = helperDotProduct(norm, lights[i].direction);
		float nDotE = helperDotProduct(norm, camera);
		/* flip normal */
		if (nDotE < 0 && nDotL < 0) {
			/* flip normal */
			GzCoord flippedNorm;
			for (int j = 0; j < 3; j++) {
				flippedNorm[j] = (-1.0f * norm[j]);
			}
			nDotL = helperDotProduct(flippedNorm, lights[i].direction);
			nDotE = helperDotProduct(flippedNorm, camera);

			/* calculate R */
			for (int j = 0; j < 3; j++) {
				R[j] = (flippedNorm[j] * 2 * nDotL) - lights[i].direction[j];
			}

			/* R dot E */
			float rDotE = helperDotProduct(R, camera);
			if (rDotE > 1) {
				rDotE = 1;
			}
			else if (rDotE < 0) {
				rDotE = 0;
			}
			float resultOfSpecPower = pow(rDotE, spec);

			for (int k = 0; k < 3; k++) {
				specColor[k] += (resultOfSpecPower * lights[i].color[k]);
				diffusedColor[k] += (nDotL * lights[i].color[k]);
			}

		}
		else if (nDotL > 0 && nDotE > 0) {

			/* calculate R */
			for (int j = 0; j < 3; j++) {
				R[j] = (norm[j] * 2 * nDotL) - lights[i].direction[j];
			}

			/* R dot E */
			float rDotE = helperDotProduct(R, camera);
			if (rDotE > 1) {
				rDotE = 1;
			}
			else if (rDotE < 0) {
				rDotE = 0;
			}
			float resultOfSpecPower = pow(rDotE, spec);
			for (int k = 0; k < 3; k++) {
				specColor[k] += (resultOfSpecPower * lights[i].color[k]);
				diffusedColor[k] += (nDotL * lights[i].color[k]);
			}
		}
	}
	/* calculate color */
	float c;
	for (int i = 0; i < 3; i++) {
		c = (ks[i] * specColor[i]) + (kd[i] * diffusedColor[i]) + (ka[i] * ambientLight.color[i]);
		if (c > 1) {
			color[i] = 1.0f;
		}
		else if (c < 0) {
			color[i] = 0;
		}
		else {
			color[i] = c;
		}
	}
}

/* Bound checking for pixel intensity to see if the color passes color limit*/
GzIntensity boundCheck(GzIntensity pixelIntensity) {
	if (pixelIntensity < 0) {
		pixelIntensity = 0;
	}
	else if (pixelIntensity > LIMIT) {
		pixelIntensity = LIMIT;
	}
	return pixelIntensity;
}

/* Bit shifting from 16 bits to 12 bits and return a character value for the color */
GzIntensity bitShiftingToChar(GzIntensity color) {
	GzIntensity shiftedColorValue = color >> 4;
	char colorChar = (char)(shiftedColorValue & 0xFF);
	return colorChar;
}


/* sorting vertices based on Y-axis */
void sortVertices(Vertex vertices[3]) {
	Vertex* temp;
	for (int i = 1; i < 3; i++) {
		Vertex keyVertex = vertices[i];
		int j = i - 1;
		while (j >= 0 && vertices[j].vertex[1] > keyVertex.vertex[1]) {
			temp = &vertices[j + 1];
			vertices[j + 1] = vertices[j];
			vertices[j] = *temp;
			j = j - 1;
		}
		vertices[j + 1] = keyVertex;
	}
}

/* To normalize vectors */
void normalize(GzCoord norm, GzCoord normalized) {
	float mag = getMagnitude(norm);
	for (int i = 0; i < 3; i++) {
		normalized[i] = norm[i] / mag;
	}
}

/* Simple implementation of scan-line algorithm with interpolations of x-y coordinates */
void drawTriangle(GzRender* r, Vertex v1, Vertex v2, Vertex v3) {

	DDA edge1_2;
	DDA edge1_3;
	DDA* left = NULL;
	DDA* right = NULL;

	float maxY = v3.vertex[1];
	float minY = v1.vertex[1];
	/* get colors of prior pixel */
	GzIntensity priorRed;
	GzIntensity priorGreen;
	GzIntensity priorBlue;
	GzIntensity priorAlpha;
	GzDepth priorZ;

	/* getting the colors */
	GzIntensity red = r->ctoi(r->flatcolor[0]);
	GzIntensity green = r->ctoi(r->flatcolor[1]);
	GzIntensity blue = r->ctoi(r->flatcolor[2]);
	GzIntensity alpha = 1;
	GzDepth z;

	/* check slope */

	edge1_2.setValues(v1, v2);
	edge1_3.setValues(v1, v3);

	/* set edges */
	if (edge1_2.slopeX < edge1_3.slopeX) {
		left = &edge1_2;
		right = &edge1_3;

	}
	else if (edge1_2.slopeX > edge1_3.slopeX) {
		left = &edge1_3;
		right = &edge1_2;
	}
	float changeInY = ceil(left->current[1]) - left->start[1];
	left->current[0] = left->start[0] + (left->slopeX * changeInY);
	right->current[0] = right->start[0] + (right->slopeX * changeInY);
	left->current[1] = left->start[1] + changeInY;
	right->current[1] = right->start[1] + changeInY;
	left->current[2] = left->start[2] + (left->slopeZ * changeInY);
	right->current[2] = right->start[2] + (right->slopeZ * changeInY);

	/* start loop */
	while (edge1_2.current[1] <= edge1_2.end[1]) {
		/* span */
		float changeInX = ceil(left->current[0]) - left->current[0];
		int startX = left->current[0] + changeInX;
		int endX = right->current[0];

		float slopeZ = (right->current[2] - left->current[2]) / (right->current[0] - left->current[0]);
		z = (GzDepth)(left->current[2] + (changeInX * slopeZ));
		int startY = edge1_2.current[1];

		while (startX <= endX) {
			/* check for z value */
			r->GzGet(startX, startY, &priorRed, &priorGreen, &priorBlue, &priorAlpha, &priorZ);
			if (priorZ > z) {
				r->GzPut(startX, startY, red, green, blue, 1, z);
			}
			startX++;
			z = z + slopeZ;
		}
		right->current[0] += right->slopeX;
		left->current[0] += left->slopeX;
		edge1_2.current[1] += 1;
		edge1_3.current[1] += 1;
		right->current[2] += right->slopeZ;
		left->current[2] += left->slopeZ;
	}
	/* swap if neccessary */
	if (edge1_2.current[1] < v3.vertex[1]) {
		edge1_2.setValues(v2, v3);
		changeInY = ceil(edge1_2.current[1]) - edge1_2.start[1];
		edge1_2.current[0] = edge1_2.start[0] + (edge1_2.slopeX * changeInY);
		edge1_2.current[1] = edge1_2.start[1] + changeInY;
		edge1_2.current[2] = edge1_2.start[2] + (edge1_2.slopeZ * changeInY);
	}

	while (edge1_2.current[1] <= edge1_2.end[1]) {
		/* span */
		float changeInX = ceil(left->current[0]) - left->current[0];
		int startX = left->current[0] + changeInX;
		int endX = right->current[0];

		float slopeZ = (right->current[2] - left->current[2]) / (right->current[0] - left->current[0]);
		z = left->current[2] + (changeInX * slopeZ);
		int startY = edge1_2.current[1];

		while (startX <= endX) {
			/* check for z value */
			r->GzGet(startX, startY, &priorRed, &priorGreen, &priorBlue, &priorAlpha, &priorZ);
			if (priorZ > z) {
				r->GzPut(startX, startY, red, green, blue, 1, z);
			}
			startX++;
			z = z + slopeZ;
		}
		right->current[0] += right->slopeX;
		left->current[0] += left->slopeX;
		edge1_2.current[1] += 1;
		edge1_3.current[1] += 1;
		right->current[2] += right->slopeZ;
		left->current[2] += left->slopeZ;
	}
}

/* Implementation of Scan-line algorithm including color interpolation */
void drawTriangleGouraudShading(GzRender* r, Vertex v1, Vertex v2, Vertex v3) {

	DDA edge1_2;
	DDA edge1_3;
	DDA* left = NULL;
	DDA* right = NULL;


	/* get colors of prior pixel */
	GzIntensity priorRed;
	GzIntensity priorGreen;
	GzIntensity priorBlue;
	GzIntensity priorAlpha;
	GzDepth priorZ;


	GzIntensity red;
	GzIntensity green;
	GzIntensity blue;
	GzDepth z;


	/* check slope */

	edge1_2.setValuesColor(v1, v2);
	edge1_3.setValuesColor(v1, v3);

	/* set edges */
	if (edge1_2.slopeX < edge1_3.slopeX) {
		left = &edge1_2;
		right = &edge1_3;
	}
	else if (edge1_2.slopeX > edge1_3.slopeX) {
		left = &edge1_3;
		right = &edge1_2;

	}
	float changeInY = ceil(left->current[1]) - left->start[1];
	left->current[0] = left->start[0] + (left->slopeX * changeInY);
	right->current[0] = right->start[0] + (right->slopeX * changeInY);
	left->current[1] = left->start[1] + changeInY;
	right->current[1] = right->start[1] + changeInY;
	left->current[2] = left->start[2] + (left->slopeZ * changeInY);
	right->current[2] = right->start[2] + (right->slopeZ * changeInY);

	/* colors */
	left->currentColor[RED] = left->startColor[RED] + (left->redSlope * changeInY);
	right->currentColor[RED] = right->startColor[RED] + (right->redSlope * changeInY);
	left->currentColor[GREEN] = left->startColor[GREEN] + (left->greenSlope * changeInY);
	right->currentColor[GREEN] = right->startColor[GREEN] + (right->greenSlope * changeInY);
	left->currentColor[BLUE] = left->startColor[BLUE] + (left->blueSlope * changeInY);
	right->currentColor[BLUE] = right->startColor[BLUE] + (right->blueSlope * changeInY);

	/* textures */
	left->currentUV[U] = left->startUV[U] + (left->slopeU * changeInY);
	right->currentUV[U] = right->startUV[U] + (right->slopeU * changeInY);
	left->currentUV[V] = left->startUV[V] + (left->slopeV * changeInY);
	right->currentUV[V] = right->startUV[V] + (right->slopeV * changeInY);

	GzColor currColor;
	GzTextureIndex currUV;
	/* start loop */
	while (edge1_2.current[1] <= edge1_2.end[1]) {
		/* span */
		float changeInX = ceil(left->current[0]) - left->current[0];
		int startX = left->current[0] + changeInX;
		int endX = right->current[0];

		float slopeZ = (right->current[2] - left->current[2]) / (right->current[0] - left->current[0]);
		/* color slopes */
		float redSlope = (right->currentColor[RED] - left->currentColor[RED]) / (right->current[0] - left->current[0]);
		float greenSlope = (right->currentColor[GREEN] - left->currentColor[GREEN]) / (right->current[0] - left->current[0]);
		float blueSlope = (right->currentColor[BLUE] - left->currentColor[BLUE]) / (right->current[0] - left->current[0]);
		z = (GzDepth)(left->current[2] + (changeInX * slopeZ));
		int startY = edge1_2.current[1];
		currColor[RED] = left->currentColor[RED];
		currColor[GREEN] = left->currentColor[GREEN];
		currColor[BLUE] = left->currentColor[BLUE];
		/* uv slopes */
		float uSlope = (right->currentUV[U] - left->currentUV[U]) / (right->current[0] - left->current[0]);
		float vSlope = (right->currentUV[V] - left->currentUV[V]) / (right->current[0] - left->current[0]);

		currUV[U] = left->currentUV[U];
		currUV[V] = left->currentUV[V];

		while (startX <= endX) {
			GzTextureIndex unwarped;

			/* check for z value */
			r->GzGet(startX, startY, &priorRed, &priorGreen, &priorBlue, &priorAlpha, &priorZ);
			if (priorZ > z) {
				unwarp(z, currUV, unwarped);
				r->tex_fun(unwarped[U], unwarped[V], currColor);
				/* getting the colors */
				red = r->ctoi(currColor[RED]);
				green = r->ctoi(currColor[GREEN]);
				blue = r->ctoi(currColor[BLUE]);
				r->GzPut(startX, startY, red, green, blue, 1, z);
			}
			startX++;
			z = z + slopeZ;
			/* color interpolation x-axis */
			currColor[RED] = currColor[RED] + redSlope;
			currColor[GREEN] = currColor[GREEN] + greenSlope;
			currColor[BLUE] = currColor[BLUE] + blueSlope;
			/* uv interpolation */
			currUV[U] += uSlope;
			currUV[V] += vSlope;

		}
		right->current[0] += right->slopeX;
		left->current[0] += left->slopeX;
		edge1_2.current[1] += 1;
		edge1_3.current[1] += 1;
		right->current[2] += right->slopeZ;
		left->current[2] += left->slopeZ;
		/* colors interpolation y-axis */
		left->currentColor[RED] += (left->redSlope);
		right->currentColor[RED] += (right->redSlope);
		left->currentColor[GREEN] += (left->greenSlope);
		right->currentColor[GREEN] += (right->greenSlope);
		left->currentColor[BLUE] += (left->blueSlope);
		right->currentColor[BLUE] += (right->blueSlope);
	}
	/* swap if neccessary */
	if (edge1_2.current[1] < v3.vertex[1]) {
		edge1_2.setValuesColor(v2, v3);
		changeInY = ceil(edge1_2.current[1]) - edge1_2.start[1];
		edge1_2.current[0] = edge1_2.start[0] + (edge1_2.slopeX * changeInY);
		edge1_2.current[1] = edge1_2.start[1] + changeInY;
		edge1_2.current[2] = edge1_2.start[2] + (edge1_2.slopeZ * changeInY);

		edge1_2.currentColor[RED] = edge1_2.startColor[RED] + (edge1_2.redSlope * changeInY);
		edge1_2.currentColor[GREEN] = edge1_2.startColor[GREEN] + (edge1_2.greenSlope * changeInY);
		edge1_2.currentColor[BLUE] = edge1_2.startColor[BLUE] + (edge1_2.blueSlope * changeInY);
	}

	while (edge1_2.current[1] <= edge1_2.end[1]) {
		/* span */
		float changeInX = ceil(left->current[0]) - left->current[0];
		int startX = left->current[0] + changeInX;
		int endX = right->current[0];

		float slopeZ = (right->current[2] - left->current[2]) / (right->current[0] - left->current[0]);
		z = left->current[2] + (changeInX * slopeZ);
		int startY = edge1_2.current[1];
		float redSlope = (right->currentColor[RED] - left->currentColor[RED]) / (right->current[0] - left->current[0]);
		float greenSlope = (right->currentColor[GREEN] - left->currentColor[GREEN]) / (right->current[0] - left->current[0]);
		float blueSlope = (right->currentColor[BLUE] - left->currentColor[BLUE]) / (right->current[0] - left->current[0]);
		currColor[RED] = left->currentColor[RED];
		currColor[GREEN] = left->currentColor[GREEN];
		currColor[BLUE] = left->currentColor[BLUE];

		while (startX <= endX) {
			/* getting the colors */
			red = r->ctoi(currColor[RED]);
			green = r->ctoi(currColor[GREEN]);
			blue = r->ctoi(currColor[BLUE]);
			/* check for z value */
			r->GzGet(startX, startY, &priorRed, &priorGreen, &priorBlue, &priorAlpha, &priorZ);
			if (priorZ > z) {
				r->GzPut(startX, startY, red, green, blue, 1, z);
			}
			startX++;
			z = z + slopeZ;
			/* color interpolation x-axis */
			currColor[RED] = currColor[RED] + redSlope;
			currColor[GREEN] = currColor[GREEN] + greenSlope;
			currColor[BLUE] = currColor[BLUE] + blueSlope;
		}
		right->current[0] += right->slopeX;
		left->current[0] += left->slopeX;
		edge1_2.current[1] += 1;
		edge1_3.current[1] += 1;
		right->current[2] += right->slopeZ;
		left->current[2] += left->slopeZ;

		/* colors interpolation y-axis */
		left->currentColor[RED] += (left->redSlope);
		right->currentColor[RED] += (right->redSlope);
		left->currentColor[GREEN] += (left->greenSlope);
		right->currentColor[GREEN] += (right->greenSlope);
		left->currentColor[BLUE] += (left->blueSlope);
		right->currentColor[BLUE] += (right->blueSlope);
	}
}

/* Implemntation of the scan line using DDA including normals interpolation */
void drawTrianglePhongShading(GzRender* r, Vertex v1, Vertex v2, Vertex v3, GzLight* lights, GzLight ambient, GzCoord camera, float spec, GzColor ks, GzColor kd, GzColor ka) {

	DDA edge1_2;
	DDA edge1_3;
	DDA* left = NULL;
	DDA* right = NULL;


	/* get colors of prior pixel */
	GzIntensity priorRed;
	GzIntensity priorGreen;
	GzIntensity priorBlue;
	GzIntensity priorAlpha;
	GzDepth priorZ;


	GzIntensity red;
	GzIntensity green;
	GzIntensity blue;
	GzDepth z;


	/* check slope */

	edge1_2.setValuesNormal(v1, v2);
	edge1_3.setValuesNormal(v1, v3);

	/* set edges */
	if (edge1_2.slopeX < edge1_3.slopeX) {
		left = &edge1_2;
		right = &edge1_3;
	}
	else if (edge1_2.slopeX > edge1_3.slopeX) {
		left = &edge1_3;
		right = &edge1_2;

	}
	float changeInY = ceil(left->current[1]) - left->start[1];
	left->current[0] = left->start[0] + (left->slopeX * changeInY);
	right->current[0] = right->start[0] + (right->slopeX * changeInY);
	left->current[1] = left->start[1] + changeInY;
	right->current[1] = right->start[1] + changeInY;
	left->current[2] = left->start[2] + (left->slopeZ * changeInY);
	right->current[2] = right->start[2] + (right->slopeZ * changeInY);

	/* normals */
	left->currentNormal[0] = left->startNormal[0] + (left->normalXSlope * changeInY);
	right->currentNormal[0] = right->startNormal[0] + (right->normalXSlope * changeInY);
	left->currentNormal[1] = left->startNormal[1] + (left->normalYSlope * changeInY);
	right->currentNormal[1] = right->startNormal[1] + (right->normalYSlope * changeInY);
	left->currentNormal[2] = left->startNormal[2] + (left->normalZSlope * changeInY);
	right->currentNormal[2] = right->startNormal[2] + (right->normalZSlope * changeInY);

	GzCoord currNormal;
	GzTextureIndex currUV;
	/* start loop */
	while (edge1_2.current[1] <= edge1_2.end[1]) {
		/* span */
		float changeInX = ceil(left->current[0]) - left->current[0];
		int startX = left->current[0] + changeInX;
		int endX = right->current[0];

		float slopeZ = (right->current[2] - left->current[2]) / (right->current[0] - left->current[0]);
		/* normal slopes */
		float xSlope = (right->currentNormal[0] - left->currentNormal[0]) / (right->current[0] - left->current[0]);
		float ySlope = (right->currentNormal[1] - left->currentNormal[1]) / (right->current[0] - left->current[0]);
		float zSlope = (right->currentNormal[2] - left->currentNormal[2]) / (right->current[0] - left->current[0]);
		
		float currZ = (left->current[2] + (changeInX * slopeZ));
		int startY = edge1_2.current[1];
		currNormal[0] = left->currentNormal[0];
		currNormal[1] = left->currentNormal[1];
		currNormal[2] = left->currentNormal[2];


		while (startX <= endX) {
			/* getting the colors */
			GzColor color;
			GzCoord normalizedNorm;

			/* normalize */
			normalize(currNormal, normalizedNorm);
			GzColor textureColor;

			/* check for z value */
			r->GzGet(startX, startY, &priorRed, &priorGreen, &priorBlue, &priorAlpha, &priorZ);
			z = (GzDepth)currZ;
			if (priorZ > z) {
				detailShade(normalizedNorm, color, lights, ambient, camera, spec, ks, kd, ka);
				red = r->ctoi(color[RED]);
				green = r->ctoi(color[GREEN]);
				blue = r->ctoi(color[BLUE]);
				r->GzPut(startX, startY, red, green, blue, 1, z);
			}
			startX++;
			currZ += slopeZ;
			/* normal interpolation x-axis */
			currNormal[0] = currNormal[0] + xSlope;
			currNormal[1] = currNormal[1] + ySlope;
			currNormal[2] = currNormal[2] + zSlope;

		}
		right->current[0] += right->slopeX;
		left->current[0] += left->slopeX;
		edge1_2.current[1] += 1;
		edge1_3.current[1] += 1;
		right->current[2] += right->slopeZ;
		left->current[2] += left->slopeZ;
		/* normals interpolation y-axis */
		left->currentNormal[0] += (left->normalXSlope);
		right->currentNormal[0] += (right->normalXSlope);
		left->currentNormal[1] += (left->normalYSlope);
		right->currentNormal[1] += (right->normalYSlope);
		left->currentNormal[2] += (left->normalZSlope);
		right->currentNormal[2] += (right->normalZSlope);

	}
	/* swap if neccessary */
	if (edge1_2.current[1] < v3.vertex[1]) {
		edge1_2.setValuesNormal(v2, v3);
		changeInY = ceil(edge1_2.current[1]) - edge1_2.start[1];
		edge1_2.current[0] = edge1_2.start[0] + (edge1_2.slopeX * changeInY);
		edge1_2.current[1] = edge1_2.start[1] + changeInY;
		edge1_2.current[2] = edge1_2.start[2] + (edge1_2.slopeZ * changeInY);

		edge1_2.currentNormal[0] = edge1_2.startNormal[0] + (edge1_2.normalXSlope * changeInY);
		edge1_2.currentNormal[1] = edge1_2.startNormal[1] + (edge1_2.normalYSlope * changeInY);
		edge1_2.currentNormal[2] = edge1_2.startNormal[2] + (edge1_2.normalZSlope * changeInY);

	}

	while (edge1_2.current[1] <= edge1_2.end[1]) {
		/* span */
		float changeInX = ceil(left->current[0]) - left->current[0];
		int startX = left->current[0] + changeInX;
		int endX = right->current[0];

		float slopeZ = (right->current[2] - left->current[2]) / (right->current[0] - left->current[0]);
		z = left->current[2] + (changeInX * slopeZ);
		int startY = edge1_2.current[1];
		/* normal slopes */
		float xSlope = (right->currentNormal[0] - left->currentNormal[0]) / (right->current[0] - left->current[0]);
		float ySlope = (right->currentNormal[1] - left->currentNormal[1]) / (right->current[0] - left->current[0]);
		float zSlope = (right->currentNormal[2] - left->currentNormal[2]) / (right->current[0] - left->current[0]);
		float currZ = left->current[2] + (changeInX * slopeZ);

		currNormal[0] = left->currentNormal[0];
		currNormal[1] = left->currentNormal[1];
		currNormal[2] = left->currentNormal[2];


		while (startX <= endX) {
			/* getting the colors */
			GzColor color;
			GzCoord normalizedNorm;

			/* normalize */
			normalize(currNormal, normalizedNorm);
			GzColor textureColor;

			/* check for z value */
			r->GzGet(startX, startY, &priorRed, &priorGreen, &priorBlue, &priorAlpha, &priorZ);
			z = (GzDepth)currZ;
			if (priorZ > z) {
				detailShade(normalizedNorm, color, lights, ambient, camera, spec, ks, kd, ka);
				red = r->ctoi(color[RED]);
				green = r->ctoi(color[GREEN]);
				blue = r->ctoi(color[BLUE]);
				r->GzPut(startX, startY, red, green, blue, 1, z);
			}
			startX++;
			currZ += slopeZ;
			/* normal interpolation x-axis */
			currNormal[0] = currNormal[0] + xSlope;
			currNormal[1] = currNormal[1] + ySlope;
			currNormal[2] = currNormal[2] + zSlope;
		}
		right->current[0] += right->slopeX;
		left->current[0] += left->slopeX;
		edge1_2.current[1] += 1;
		edge1_3.current[1] += 1;
		right->current[2] += right->slopeZ;
		left->current[2] += left->slopeZ;

		/* normal interpolation y-axis */
		left->currentNormal[0] += (left->normalXSlope);
		right->currentNormal[0] += (right->normalXSlope);
		left->currentNormal[1] += (left->normalYSlope);
		right->currentNormal[1] += (right->normalYSlope);
		left->currentNormal[2] += (left->normalZSlope);
		right->currentNormal[2] += (right->normalZSlope);
	}
}
/* mode 1 -> Phong, mode 2 -> Gouraud 
	This function helps in shading with a texture map */
void drawTriangleWithTexture(int mode, GzRender* r, Vertex v1, Vertex v2, Vertex v3, GzLight* lights, GzLight ambient, GzCoord camera, float spec, GzColor ks) {

	DDA edge1_2;
	DDA edge1_3;
	DDA* left = NULL;
	DDA* right = NULL;


	/* get colors of prior pixel */
	GzIntensity priorRed;
	GzIntensity priorGreen;
	GzIntensity priorBlue;
	GzIntensity priorAlpha;
	GzDepth priorZ;


	GzIntensity red;
	GzIntensity green;
	GzIntensity blue;
	GzDepth z;


	/* check slope */

	edge1_2.setValuesNormal(v1, v2);
	edge1_3.setValuesNormal(v1, v3);

	/* set edges */
	if (edge1_2.slopeX < edge1_3.slopeX) {
		left = &edge1_2;
		right = &edge1_3;
	}
	else if (edge1_2.slopeX > edge1_3.slopeX) {
		left = &edge1_3;
		right = &edge1_2;

	}
	float changeInY = ceil(left->current[1]) - left->start[1];
	left->current[0] = left->start[0] + (left->slopeX * changeInY);
	right->current[0] = right->start[0] + (right->slopeX * changeInY);
	left->current[1] = left->start[1] + changeInY;
	right->current[1] = right->start[1] + changeInY;
	left->current[2] = left->start[2] + (left->slopeZ * changeInY);
	right->current[2] = right->start[2] + (right->slopeZ * changeInY);

	/* normals */
	left->currentNormal[0] = left->startNormal[0] + (left->normalXSlope * changeInY);
	right->currentNormal[0] = right->startNormal[0] + (right->normalXSlope * changeInY);
	left->currentNormal[1] = left->startNormal[1] + (left->normalYSlope * changeInY);
	right->currentNormal[1] = right->startNormal[1] + (right->normalYSlope * changeInY);
	left->currentNormal[2] = left->startNormal[2] + (left->normalZSlope * changeInY);
	right->currentNormal[2] = right->startNormal[2] + (right->normalZSlope * changeInY);

	/* textures */
	//left->currentUV[U] = left->startUV[U] + (left->slopeU * changeInY);
	//right->currentUV[U] = right->startUV[U] + (right->slopeU * changeInY);
	//left->currentUV[V] = left->startUV[V] + (left->slopeV * changeInY);
	//right->currentUV[V] = right->startUV[V] + (right->slopeV * changeInY);

	GzCoord currNormal;
	GzTextureIndex currUV;
	/* start loop */
	while (edge1_2.current[1] <= edge1_2.end[1]) {
		/* span */
		float changeInX = ceil(left->current[0]) - left->current[0];
		int startX = left->current[0] + changeInX;
		int endX = right->current[0];

		float slopeZ = (right->current[2] - left->current[2]) / (right->current[0] - left->current[0]);
		/* normal slopes */
		float xSlope = (right->currentNormal[0] - left->currentNormal[0]) / (right->current[0] - left->current[0]);
		float ySlope = (right->currentNormal[1] - left->currentNormal[1]) / (right->current[0] - left->current[0]);
		float zSlope = (right->currentNormal[2] - left->currentNormal[2]) / (right->current[0] - left->current[0]);

		float currZ = (left->current[2] + (changeInX * slopeZ));
		int startY = edge1_2.current[1];
		currNormal[0] = left->currentNormal[0];
		currNormal[1] = left->currentNormal[1];
		currNormal[2] = left->currentNormal[2];

		currUV[U] = left->currentUV[U];
		currUV[V] = left->currentUV[V];

		/* uv slopes */
		float uSlope = (right->currentUV[U] - left->currentUV[U]) / (right->current[0] - left->current[0]);
		float vSlope = (right->currentUV[V] - left->currentUV[V]) / (right->current[0] - left->current[0]);

		currUV[U] = left->currentUV[U] + (uSlope * changeInX);
		currUV[V] = left->currentUV[V] + (vSlope * changeInX);


		while (startX <= endX) {
			/* getting the colors */
			GzColor color;
			GzCoord normalizedNorm;
			GzTextureIndex unwarpedUV;

			/* normalize */
			normalize(currNormal, normalizedNorm);
			GzColor textureColor;

			/* check for z value */
			r->GzGet(startX, startY, &priorRed, &priorGreen, &priorBlue, &priorAlpha, &priorZ);
			z = (GzDepth)currZ;
			if (priorZ > z) {
				unwarp(currZ, currUV, unwarpedUV);
				r->tex_fun(unwarpedUV[U], unwarpedUV[V], textureColor);
				if (mode == 1) {
					textureShadePhong(normalizedNorm, color, lights, camera, ks, textureColor, ambient, spec);
				}
				else if (mode == 2) {
					textureShadeGouraud(normalizedNorm, color, lights, camera, textureColor, ambient, spec);
				}

				red = r->ctoi(color[RED]);
				green = r->ctoi(color[GREEN]);
				blue = r->ctoi(color[BLUE]);
				r->GzPut(startX, startY, red, green, blue, 1, z);
			}
			startX++;
			currZ += slopeZ;
			/* normal interpolation x-axis */
			currNormal[0] = currNormal[0] + xSlope;
			currNormal[1] = currNormal[1] + ySlope;
			currNormal[2] = currNormal[2] + zSlope;
			/* uv interpolation */
			currUV[U] += uSlope;
			currUV[V] += vSlope;

		}
		right->current[0] += right->slopeX;
		left->current[0] += left->slopeX;
		edge1_2.current[1] += 1;
		edge1_3.current[1] += 1;
		right->current[2] += right->slopeZ;
		left->current[2] += left->slopeZ;
		/* normals interpolation y-axis */
		left->currentNormal[0] += (left->normalXSlope);
		right->currentNormal[0] += (right->normalXSlope);
		left->currentNormal[1] += (left->normalYSlope);
		right->currentNormal[1] += (right->normalYSlope);
		left->currentNormal[2] += (left->normalZSlope);
		right->currentNormal[2] += (right->normalZSlope);
		/* uv interpolation y-axis */
		left->currentUV[U] += (left->slopeU);
		right->currentUV[U] += (right->slopeU);
		left->currentUV[V] += (left->slopeV);
		right->currentUV[V] += (right->slopeV);
	}
	/* swap if neccessary */
	if (edge1_2.current[1] < v3.vertex[1]) {
		edge1_2.setValuesNormal(v2, v3);
		changeInY = ceil(edge1_2.current[1]) - edge1_2.start[1];
		edge1_2.current[0] = edge1_2.start[0] + (edge1_2.slopeX * changeInY);
		edge1_2.current[1] = edge1_2.start[1] + changeInY;
		edge1_2.current[2] = edge1_2.start[2] + (edge1_2.slopeZ * changeInY);

		edge1_2.currentNormal[0] = edge1_2.startNormal[0] + (edge1_2.normalXSlope * changeInY);
		edge1_2.currentNormal[1] = edge1_2.startNormal[1] + (edge1_2.normalYSlope * changeInY);
		edge1_2.currentNormal[2] = edge1_2.startNormal[2] + (edge1_2.normalZSlope * changeInY);

		//edge1_2.currentUV[U] = edge1_2.startUV[U] + (edge1_2.slopeU * changeInY);
		//edge1_2.currentUV[V] = edge1_2.startUV[V] + (edge1_2.slopeV * changeInY);
	}

	while (edge1_2.current[1] <= edge1_2.end[1]) {
		/* span */
		float changeInX = ceil(left->current[0]) - left->current[0];
		int startX = left->current[0] + changeInX;
		int endX = right->current[0];

		float slopeZ = (right->current[2] - left->current[2]) / (right->current[0] - left->current[0]);
		z = left->current[2] + (changeInX * slopeZ);
		int startY = edge1_2.current[1];
		/* normal slopes */
		float xSlope = (right->currentNormal[0] - left->currentNormal[0]) / (right->current[0] - left->current[0]);
		float ySlope = (right->currentNormal[1] - left->currentNormal[1]) / (right->current[0] - left->current[0]);
		float zSlope = (right->currentNormal[2] - left->currentNormal[2]) / (right->current[0] - left->current[0]);
		float currZ = left->current[2] + (changeInX * slopeZ);

		currNormal[0] = left->currentNormal[0];
		currNormal[1] = left->currentNormal[1];
		currNormal[2] = left->currentNormal[2];

		/* uv slopes */
		float uSlope = (right->currentUV[U] - left->currentUV[U]) / (right->current[0] - left->current[0]);
		float vSlope = (right->currentUV[V] - left->currentUV[V]) / (right->current[0] - left->current[0]);


		currUV[U] = left->currentUV[U] + (uSlope * changeInX);
		currUV[V] = left->currentUV[V] + (vSlope * changeInX);

		while (startX <= endX) {
			/* getting the colors */
			GzColor color;
			GzCoord normalizedNorm;
			GzTextureIndex unwarpedUV;

			/* normalize */
			normalize(currNormal, normalizedNorm);
			GzColor textureColor;

			/* check for z value */
			r->GzGet(startX, startY, &priorRed, &priorGreen, &priorBlue, &priorAlpha, &priorZ);
			z = (GzDepth)currZ;
			if (priorZ > z) {
				unwarp(currZ, currUV, unwarpedUV);
				r->tex_fun(unwarpedUV[U], unwarpedUV[V], textureColor);
				if (mode == 1) {
					textureShadePhong(normalizedNorm, color, lights, camera, ks, textureColor, ambient, spec);
				}
				else if (mode == 2) {
					textureShadeGouraud(normalizedNorm, color, lights, camera, textureColor, ambient, spec);
				}
				red = r->ctoi(color[RED]);
				green = r->ctoi(color[GREEN]);
				blue = r->ctoi(color[BLUE]);
				r->GzPut(startX, startY, red, green, blue, 1, z);
			}
			startX++;
			currZ += slopeZ;
			/* normal interpolation x-axis */
			currNormal[0] = currNormal[0] + xSlope;
			currNormal[1] = currNormal[1] + ySlope;
			currNormal[2] = currNormal[2] + zSlope;
			/* uv interpolation */
			currUV[U] += uSlope;
			currUV[V] += vSlope;

		}
		right->current[0] += right->slopeX;
		left->current[0] += left->slopeX;
		edge1_2.current[1] += 1;
		edge1_3.current[1] += 1;
		right->current[2] += right->slopeZ;
		left->current[2] += left->slopeZ;

		/* normal interpolation y-axis */
		left->currentNormal[0] += (left->normalXSlope);
		right->currentNormal[0] += (right->normalXSlope);
		left->currentNormal[1] += (left->normalYSlope);
		right->currentNormal[1] += (right->normalYSlope);
		left->currentNormal[2] += (left->normalZSlope);
		right->currentNormal[2] += (right->normalZSlope);
		/* uv interpolation y-axis */
		left->currentUV[U] += (left->slopeU);
		right->currentUV[U] += (right->slopeU);
		left->currentUV[V] += (left->slopeV);
		right->currentUV[V] += (right->slopeV);
	}
}
/*warping uv coordinates */
void warp(Vertex* vert) {
	float z = vert->vertex[2];
	float vz = (float)z / ((float)MAXINT - (float)z);
	vert->texture[U] = vert->texture[U] / (vz + 1.0f);
	vert->texture[V] = vert->texture[V] / (vz + 1.0f);
}

/*unwarping uv coordinates */
void unwarp(float z, GzTextureIndex warped, GzTextureIndex unwarped) {
	float vz =  (float)z / ((float)MAXINT - (float)z);
	unwarped[U] = warped[U] * (vz + 1.0f);
	unwarped[V] = warped[V] * (vz + 1.0f);
}
/* calculate cosine in radians  */
float cosine(float degree) {
	return cos((double)degree * (PI / 180.0));
}

/* calculates sine in radians */
float sine(float degree) {
	return sin((double)degree * (PI / 180.0));
}

/* Does the dot-product with TOS (Top of Stack) */
void dotProduct(GzMatrix matVert, GzMatrix transformationStack[], int level) {
	for (int i = 0; i < 4; i++) {
		for (int k = 0; k < 4; k++) {
			float sum = 0.0;
			for (int j = 0; j < 4; j++) {
				sum += transformationStack[level - 1][k][j] * matVert[j][i];
			}
			transformationStack[level][k][i] = sum;
		}
	}
}

/* Gets the dissimarity between vectors (to get the orthogonal vector) */
void getCrossProduct(float* v1, float* v2, float* result) {
	result[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
	result[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
	result[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
}

/* makes pure (unitary matrix) preprocessing matrix for norms-stack */
void makePureRotationMatrix(GzMatrix regMatrix, GzMatrix pureMatrix) {
	float sampleRow[3];
	for (int i = 0; i < 3; i++) {
		sampleRow[i] = regMatrix[1][i];
	}
	float scaleFactor = getMagnitude(sampleRow);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (j == 3 && i == 3) {
				pureMatrix[i][j] = 1;
			}
			else if (i < 3 && j < 3) {
				pureMatrix[i][j] = regMatrix[i][j] / scaleFactor;
			}
			else {
				pureMatrix[i][j] = 0;
			}
		}
	}
}


/* Render Library implementation  */
GzRender::~GzRender()
{
	/* HW1.2 clean up, free buffer memory */
	delete[] framebuffer;
	delete[] pixelbuffer;

}

int GzRender::GzDefault()
{
	/* HW1.3 set pixel buffer to some default values - start a new frame */
	GzPixel defaultPixel = { 4000, 4000, 4000, 1, MAXINT };
	int resolution = xres * yres;
	for (int i = 0; i < resolution; i++) {
		pixelbuffer[i] = defaultPixel;
	}
	return GZ_SUCCESS;
}


int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	/* HW1.4 write pixel values into the buffer */
	if (i >= 0 && i < xres && j >= 0 && j < yres) {
		int index = ARRAY(i, j);
		GzPixel olderPixel = pixelbuffer[index];
		/* check if z is closer, if so, re-paint the pixel.*/
		if (olderPixel.z > z) {
			GzPixel pixel = { r, g, b, a, z };
			pixelbuffer[index] = pixel;
		}
	}

	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z)
{
	/* HW1.5 retrieve a pixel information from the pixel buffer */
	if (i >= 0 && i < xres && j >= 0 && j < yres) {
		int index = ARRAY(i, j);
		GzPixel pixel = pixelbuffer[index];
		*r = pixel.red;
		*g = pixel.green;
		*b = pixel.blue;
		*a = pixel.alpha;
		*z = pixel.z;
	}

	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
	/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */
	fprintf(outfile, "P6 %d %d 255\n", xres, yres);
	int resolution = xres * yres;
	for (int i = 0; i < resolution; i++) {
		GzPixel pixel = pixelbuffer[i];
		GzIntensity red = boundCheck(pixel.red);
		GzIntensity green = boundCheck(pixel.green);
		GzIntensity blue = boundCheck(pixel.blue);

		char r = bitShiftingToChar(red);
		char g = bitShiftingToChar(green);
		char b = bitShiftingToChar(blue);

		char color[3];
		color[0] = r;
		color[1] = g;
		color[2] = b;
		fwrite(color, 1, 3, outfile);
	}

	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
	/* HW1.7 write pixels to framebuffer:
		- put the pixels into the frame buffer
		- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
		- NOT red, green, and blue !!!
	*/
	int resolution = xres * yres;
	for (int i = 0; i < resolution; i++) {
		GzPixel pixel = pixelbuffer[i];
		GzIntensity red = boundCheck(pixel.red);
		GzIntensity green = boundCheck(pixel.green);
		GzIntensity blue = boundCheck(pixel.blue);

		char r = bitShiftingToChar(red);
		char g = bitShiftingToChar(green);
		char b = bitShiftingToChar(blue);

		framebuffer[RGB * i] = b;
		framebuffer[RGB * i + 1] = g;
		framebuffer[RGB * i + 2] = r;
	}

	return GZ_SUCCESS;
}




/***********************************************/
/* HW2 methods: implement from here */

int GzRender::GzPutAttribute(int numAttributes, GzToken* nameList, GzPointer* valueList)
{
	/* HW 2.1
	-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	-- In later homeworks set shaders, interpolaters, texture maps, and lights
	*/

	if (numAttributes == 2 && nameList[0] == GZ_AASHIFTX) {
		float* xShift = (float*)valueList[0];
		float* yShift = (float*)valueList[1];
		/* getting the offsets */
		offset[0] = *xShift;
		offset[1] = *yShift;
	}
	else if (numAttributes == 1 && nameList[0] == GZ_AMBIENT_LIGHT) {
		GzLight* ambient = (GzLight*)valueList[0];
		ambientlight.direction[0] = ambient->direction[0];
		ambientlight.direction[1] = ambient->direction[1];
		ambientlight.direction[2] = ambient->direction[2];

		ambientlight.color[RED] = ambient->color[RED];
		ambientlight.color[GREEN] = ambient->color[GREEN];
		ambientlight.color[BLUE] = ambient->color[BLUE];
	}
	else if (numAttributes == 3) {
		GzLight* lightsFromApp1 = (GzLight*)valueList[0];
		GzLight* lightsFromApp2 = (GzLight*)valueList[1];
		GzLight* lightsFromApp3 = (GzLight*)valueList[2];

		lights[0].color[RED] = lightsFromApp1->color[RED];
		lights[0].color[GREEN] = lightsFromApp1->color[GREEN];
		lights[0].color[BLUE] = lightsFromApp1->color[BLUE];

		lights[0].direction[0] = lightsFromApp1->direction[0];
		lights[0].direction[1] = lightsFromApp1->direction[1];
		lights[0].direction[2] = lightsFromApp1->direction[2];



		lights[1].color[RED] = lightsFromApp2->color[RED];
		lights[1].color[GREEN] = lightsFromApp2->color[GREEN];
		lights[1].color[BLUE] = lightsFromApp2->color[BLUE];

		lights[1].direction[0] = lightsFromApp2->direction[0];
		lights[1].direction[1] = lightsFromApp2->direction[1];
		lights[1].direction[2] = lightsFromApp2->direction[2];


		lights[2].color[RED] = lightsFromApp3->color[RED];
		lights[2].color[GREEN] = lightsFromApp3->color[GREEN];
		lights[2].color[BLUE] = lightsFromApp3->color[BLUE];

		lights[2].direction[0] = lightsFromApp3->direction[0];
		lights[2].direction[1] = lightsFromApp3->direction[1];
		lights[2].direction[2] = lightsFromApp3->direction[2];
	}
	else if (numAttributes == 6) {
		float* kd = (float*)valueList[0];
		int* style = (int*)valueList[1];
		float* ka = (float*)valueList[2];
		float* ks = (float*)valueList[3];
		float* specular = (float*)valueList[4];
		tex_fun = (GzTexture)valueList[5];

		spec = *specular;

		Kd[RED] = kd[RED];
		Kd[GREEN] = kd[GREEN];
		Kd[BLUE] = kd[BLUE];

		Ka[RED] = ka[RED];
		Ka[GREEN] = ka[GREEN];
		Ka[BLUE] = ka[BLUE];

		Ks[RED] = ks[RED];
		Ks[GREEN] = ks[GREEN];
		Ks[BLUE] = ks[BLUE];

		interp_mode = *style;
	}

	return GZ_SUCCESS;
}

int GzRender::GzPutTriangle(int	numParts, GzToken* nameList, GzPointer* valueList)
/* numParts - how many names and values */
{
	/* HW 2.2
	-- Pass in a triangle description with tokens and values corresponding to
		  GZ_NULL_TOKEN:		do nothing - no values
		  GZ_POSITION:		3 vert positions in model space
	-- Invoke the rastrizer/scanline framework
	-- Return error code
	*/
	GzTextureIndex* uvList;
	if (nameList[2] == GZ_TEXTURE_INDEX) {
		uvList = (GzTextureIndex*)valueList[2];

	}

	if (interp_mode == GZ_FLAT) {
		GzCoord E = { 0.0, 0.0, -1.0f };

		GzCoord* norm = (GzCoord*)valueList[1];
		GzCoord* pointer = (GzCoord*)valueList[0];
		/* apply the transformation */
		/* convert to homogenous matrix and apply transfromation */
		float transformedVertices4D[4][4], transformedVertices3D[3][3], transformedNormals4D[4][4], transformedNormals3D[3][3];
		for (int i = 0; i < 4; i++) {
			for (int k = 0; k < 4; k++) {
				float sum = 0.0;
				float normSum = 0.0;
				for (int j = 0; j < 4; j++) {
					if (j < 3) {
						sum += Ximage[matlevel - 1][k][j] * pointer[i][j];
						normSum += Xnorm[normlevel - 1][k][j] * norm[i][j];
					}
					else {
						sum += Ximage[matlevel - 1][k][j] * 1;
						normSum += Xnorm[normlevel - 1][k][j] * 1;
					}
				}
				transformedVertices4D[i][k] = sum;
				transformedNormals4D[i][k] = normSum;
			}
		}
		/* covert back to 3D */
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				transformedVertices3D[i][j] = (float)(transformedVertices4D[i][j] / transformedVertices4D[i][3]);
				transformedNormals3D[i][j] = (float)(transformedNormals4D[i][j] / transformedNormals4D[i][3]);

			}
		}

		/* applying the offsets */
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 2; j++) {
				transformedVertices3D[i][j] = transformedVertices3D[i][j] - offset[j];
			}
		}
		/* point back to a 3D matrix */
		pointer = transformedVertices3D;
		/* transfomed normals back to 3D */
		norm = transformedNormals3D;
		GzToken token = nameList[0];
		detailShade(norm[0], flatcolor, lights, ambientlight, E, spec, Ks, Kd, Ka);

		/* sort the vertices with respect to y-axis */
		Vertex* vert1 = new Vertex();
		Vertex* vert2 = new Vertex();
		Vertex* vert3 = new Vertex();

		float v1[3] = { pointer[0][0], pointer[0][1], pointer[0][2] };
		float v2[3] = { pointer[1][0], pointer[1][1], pointer[1][2] };
		float v3[3] = { pointer[2][0], pointer[2][1], pointer[2][2] };

		vert1->setValues(v1);
		vert2->setValues(v2);
		vert3->setValues(v3);



		Vertex container[3] = { *v1, *v2, *v3 };
		sortVertices(container);

		drawTriangle(this, container[0], container[1], container[2]);



		return GZ_SUCCESS;

	}
	else if (interp_mode == GZ_COLOR) {

		GzColor color1, color2, color3;
		GzCoord* norm = (GzCoord*)valueList[1];
		GzCoord E = { 0.0, 0.0, -1.0f };

		GzCoord* pointer = (GzCoord*)valueList[0];

		/* apply the transformation */
		/* convert to homogenous matrix and apply transfromation */
		float transformedVertices4D[4][4], transformedVertices3D[3][3], transformedNormals4D[4][4], transformedNormals3D[3][3];
		for (int i = 0; i < 4; i++) {
			for (int k = 0; k < 4; k++) {
				float sum = 0.0;
				float normSum = 0.0;
				for (int j = 0; j < 4; j++) {
					if (j < 3) {
						sum += Ximage[matlevel - 1][k][j] * pointer[i][j];
						normSum += Xnorm[normlevel - 1][k][j] * norm[i][j];
					}
					else {
						sum += Ximage[matlevel - 1][k][j] * 1;
						normSum += Xnorm[normlevel - 1][k][j] * 1;
					}
				}
				transformedVertices4D[i][k] = sum;
				transformedNormals4D[i][k] = normSum;
			}
		}
		/* covert back to 3D */
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				transformedVertices3D[i][j] = (float)(transformedVertices4D[i][j] / transformedVertices4D[i][3]);
				transformedNormals3D[i][j] = (float)(transformedNormals4D[i][j] / transformedNormals4D[i][3]);

			}
		}
		/* point back to a 3D matrix */
		pointer = transformedVertices3D;
		/* transfomed normals back to 3D */
		norm = transformedNormals3D;
		GzToken token = nameList[0];

		GzColor textureColor1;
		GzColor textureColor2;
		GzColor textureColor3;
		GzTextureIndex unwarped;


		
		/* to Gouraud shade without texture */
		//detailShade(norm[0], color1, lights, ambientlight, E, spec, Ks, Kd, Ka);
		//detailShade(norm[1], color2, lights, ambientlight, E, spec, Ks, Kd, Ka);
		//detailShade(norm[2], color3, lights, ambientlight, E, spec, Ks, Kd, Ka);

		/* sort the vertices with respect to y-axis */

		float v1[3] = { pointer[0][0], pointer[0][1], pointer[0][2] };
		float v2[3] = { pointer[1][0], pointer[1][1], pointer[1][2] };
		float v3[3] = { pointer[2][0], pointer[2][1], pointer[2][2] };

		Vertex* vert1 = new Vertex();
		Vertex* vert2 = new Vertex();
		Vertex* vert3 = new Vertex();

		GzTextureIndex uv1 = { uvList[0][0], uvList[0][1] };
		GzTextureIndex uv2 = { uvList[1][0], uvList[1][1] };
		GzTextureIndex uv3 = { uvList[2][0], uvList[2][1] };

		vert1->setValuesNormal(norm[0], v1, uv1);
		vert2->setValuesNormal(norm[1], v2, uv2);
		vert3->setValuesNormal(norm[2], v3, uv3);

		warp(vert1);
		warp(vert2);
		warp(vert3);

		Vertex container[3] = { *vert1, *vert2, *vert3 };
		sortVertices(container);
		GzCoord camera = { 0.0, 0.0, -1.0f };

		/* mode 2 for Gouraud */
		drawTriangleWithTexture(2,  this, container[0], container[1], container[2], lights, ambientlight, camera, spec, Ks);

		return GZ_SUCCESS;

	}
	else {
		GzCoord* norm = (GzCoord*)valueList[1];
		GzCoord* pointer = (GzCoord*)valueList[0];
		/* apply the transformation */
		/* convert to homogenous matrix and apply transfromation */
		float transformedVertices4D[4][4], transformedVertices3D[3][3], transformedNormals4D[4][4], transformedNormals3D[3][3];
		for (int i = 0; i < 4; i++) {
			for (int k = 0; k < 4; k++) {
				float sum = 0.0;
				float normSum = 0.0;
				for (int j = 0; j < 4; j++) {
					if (j < 3) {
						sum += Ximage[matlevel - 1][k][j] * pointer[i][j];
						normSum += Xnorm[normlevel - 1][k][j] * norm[i][j];
					}
					else {
						sum += Ximage[matlevel - 1][k][j] * 1;
						normSum += Xnorm[normlevel - 1][k][j] * 1;
					}
				}
				transformedVertices4D[i][k] = sum;
				transformedNormals4D[i][k] = normSum;
			}
		}
		/* covert back to 3D */
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				transformedVertices3D[i][j] = (float)(transformedVertices4D[i][j] / transformedVertices4D[i][3]);
				transformedNormals3D[i][j] = (float)(transformedNormals4D[i][j] / transformedNormals4D[i][3]);

			}
		}

		/* applying the offsets */
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 2; j++) {
				transformedVertices3D[i][j] = transformedVertices3D[i][j] - offset[j];
			}
		}
		/* point back to a 3D matrix */
		pointer = transformedVertices3D;
		/* transfomed normals back to 3D */
		norm = transformedNormals3D;
		GzToken token = nameList[0];

		/* sort the vertices with respect to y-axis */

		float* v1 = (float*)pointer[0];
		float* v2 = (float*)pointer[1];
		float* v3 = (float*)pointer[2];

		
		GzTextureIndex uv1 = { uvList[0][0], uvList[0][1] };
		GzTextureIndex uv2 = { uvList[1][0], uvList[1][1] };
		GzTextureIndex uv3 = { uvList[2][0], uvList[2][1] };


		Vertex* vert1 = new Vertex();
		Vertex * vert2 = new Vertex();
		Vertex* vert3 = new Vertex();

		vert1->setValuesNormal(norm[0], v1, uv1);
		vert2->setValuesNormal(norm[1], v2, uv2);
		vert3->setValuesNormal(norm[2], v3, uv3);

		warp(vert1);
		warp(vert2);
		warp(vert3);

		Vertex container[3] = { *vert1, *vert2, *vert3 };
		sortVertices(container);
		GzCoord camera = { 0.0, 0.0, -1.0f };
		/* mode 1 for Phong */
		drawTriangleWithTexture(1, this, container[0], container[1], container[2], lights, ambientlight, camera, spec, Ks);

	}


	return GZ_SUCCESS;
}



int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
	/* HW 3.1
	// Create rotate matrix : rotate along x axis
	// Pass back the matrix using mat value
	*/
	/* sine and cosines */
	float cosValue = cosine(degree);
	float sineValue = sine(degree);

	/* Building the matrices */
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				mat[i][j] = 1.0;
			}
			else {
				mat[i][j] = 0.0;
			}

		}
	}
	mat[0][0] = 1;
	mat[1][1] = cosValue;
	mat[1][2] = -sineValue;
	mat[2][1] = sineValue;
	mat[2][2] = cosValue;

	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
	/* HW 3.2
	// Create rotate matrix : rotate along y axis
	// Pass back the matrix using mat value
	*/
	/* sine and cosines */
	float cosValue = cosine(degree);
	float sineValue = sine(degree);

	/* Building the matrices */
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				mat[i][j] = 1.0;
			}
			else {
				mat[i][j] = 0.0;
			}
		}
	}
	mat[0][0] = cosValue;
	mat[1][1] = 1;
	mat[0][2] = sineValue;
	mat[2][0] = -sineValue;
	mat[2][2] = cosValue;

	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
	/* HW 3.3
	// Create rotate matrix : rotate along z axis
	// Pass back the matrix using mat value
	*/
	/* sine and cosines */
	float cosValue = cosine(degree);
	float sineValue = sine(degree);

	/* Building the matrices */
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				mat[i][j] = 1.0;
			}
			else {
				mat[i][j] = 0.0;
			}
		}
	}
	mat[0][0] = cosValue;
	mat[0][1] = -sineValue;
	mat[1][0] = sineValue;
	mat[1][1] = cosValue;
	mat[2][2] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
	/* HW 3.4
	// Create translation matrix
	// Pass back the matrix using mat value
	*/
	/* Building Matrix */
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				mat[i][j] = 1.0;
			}
			else {
				mat[i][j] = 0.0;
			}
		}
	}
	mat[0][3] = translate[0];
	mat[1][3] = translate[1];
	mat[2][3] = translate[2];


	return GZ_SUCCESS;
}


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
	/* HW 3.5
	// Create scaling matrix
	// Pass back the matrix using mat value
	*/
	/* Building matrix */
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				mat[i][j] = 1.0;
			}
			else {
				mat[i][j] = 0.0;
			}
		}
	}
	mat[0][0] = scale[0];
	mat[1][1] = scale[1];
	mat[2][2] = scale[2];

	return GZ_SUCCESS;
}


GzRender::GzRender(int xRes, int yRes)
{
	/* HW1.1 create a framebuffer for MS Windows display:
	 -- set display resolution
	 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
	 -- allocate memory for pixel buffer
	 */
	framebuffer = (char*)malloc(3 * sizeof(char) * xRes * yRes);
	xres = (unsigned short)xRes;
	yres = (unsigned short)yRes;

	int frameResolution = xres * yres;
	int frameBufferDepth = RGB * frameResolution;
	//framebuffer = new char[frameBufferDepth];
	pixelbuffer = new GzPixel[frameResolution];


	/* HW 3.6
	- setup Xsp and anything only done once
	- init default camera
	*/
	/* setting up Xsp*/
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				Xsp[i][j] = 1.0;
			}
			else {
				Xsp[i][j] = 0.0;
			}
		}
	}
	Xsp[0][0] = (float)xres / 2;
	Xsp[1][1] = (float)-1.0f * (yres / 2);
	Xsp[0][3] = (float)xres / 2;
	Xsp[1][3] = (float)yres / 2;
	Xsp[2][2] = (float)MAXINT;

	/* setting up the default camera */
	/* field of View */
	m_camera.FOV = DEFAULT_FOV;

	/* default look - at coordinates */
	m_camera.lookat[0] = 0;
	m_camera.lookat[1] = 0;
	m_camera.lookat[2] = 0;

	/* default worldup coordinates */
	m_camera.worldup[0] = 0;
	m_camera.worldup[1] = 1;
	m_camera.worldup[2] = 0;

	/* default image - plane position of camera */
	m_camera.position[0] = DEFAULT_IM_X;
	m_camera.position[1] = DEFAULT_IM_Y;
	m_camera.position[2] = DEFAULT_IM_Z;
}


int GzRender::GzBeginRender()
{
	/* HW 3.7
	- setup for start of each frame - init frame buffer color,alpha,z
	- compute Xiw and projection xform Xpi from camera definition
	- init Ximage - put Xsp at base of stack, push on Xpi and Xiw
	- now stack contains Xsw and app can push model Xforms when needed
	*/

	/* compute Xiw and Xpi */
	/* Xiw */
	GzCoord cameraX, cameraY, cameraZ;


	/* setting z-axis */
	for (int i = 0; i < 3; i++) {
		cameraZ[i] = m_camera.lookat[i] - m_camera.position[i];
	}
	float lookAtMagnitude = getMagnitude(cameraZ);

	for (int i = 0; i < 3; i++) {
		cameraZ[i] = cameraZ[i] / lookAtMagnitude;
	}


	/* setting y-axis */
	float upDiff = helperDotProduct(m_camera.worldup, cameraZ);
	for (int i = 0; i < 3; i++) {
		cameraY[i] = m_camera.worldup[i] - (upDiff * cameraZ[i]);
	}
	float realYMagnitude = getMagnitude(cameraY);

	for (int i = 0; i < 3; i++) {
		cameraY[i] = cameraY[i] / realYMagnitude;
	}


	/* setting x-axis */
	getCrossProduct(cameraY, cameraZ, cameraX);

	/* building Xiw matrix */
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == 0) {
				if (j < 3) {
					m_camera.Xiw[i][j] = cameraX[j];
				}
			}
			else if (i == 1) {
				if (j < 3) {
					m_camera.Xiw[i][j] = cameraY[j];
				}
			}
			else if (i == 2) {
				if (j < 3) {
					m_camera.Xiw[i][j] = cameraZ[j];
				}
			}
			else {
				if (j == 3) {
					m_camera.Xiw[i][j] = 1;
				}
				else {
					m_camera.Xiw[i][j] = 0;
				}
			}

		}
	}

	m_camera.Xiw[0][3] = -1.0f * helperDotProduct(cameraX, m_camera.position);
	m_camera.Xiw[1][3] = -1.0f * helperDotProduct(cameraY, m_camera.position);
	m_camera.Xiw[2][3] = -1.0f * helperDotProduct(cameraZ, m_camera.position);


	/* Xpi */
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				m_camera.Xpi[i][j] = 1;
			}
			else {
				m_camera.Xpi[i][j] = 0;
			}
		}
	}
	float focalDepth = max(0, tan((m_camera.FOV * (PI / 180.0) / 2)));
	m_camera.Xpi[2][2] = focalDepth;
	m_camera.Xpi[3][2] = focalDepth;

	/* setting up transformation stack with Xsp Xpi and Xiw */
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Ximage[0][i][j] = Xsp[i][j];
		}
	}
	int status = 0;
	matlevel = 1;

	/* setting up normalstack stack with identity matrix */
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				Xnorm[0][i][j] = 1;
			}
			else {
				Xnorm[0][i][j] = 0;
			}
		}

	}
	normlevel = 1;

	status = GzPushMatrix(m_camera.Xpi);
	if (status == 1) {
		return GZ_FAILURE;
	}
	status = GzPushMatrix(m_camera.Xiw);
	if (status == 1) {
		return GZ_FAILURE;
	}
	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
	/* HW 3.8
	/*- overwrite renderer camera structure with new camera definition
	*/
	/* updating the field of view */
	m_camera.FOV = camera.FOV;

	/* updating the postion, lookat (z) and wordup (y) vectors */
	for (int i = 0; i <= 2; i++) {
		m_camera.position[i] = camera.position[i];
		m_camera.lookat[i] = camera.lookat[i];
		m_camera.worldup[i] = camera.worldup[i];
	}

	/* updating the image - perspective and word - image transfromation matrices */
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			m_camera.Xiw[i][j] = camera.Xiw[i][j];
			m_camera.Xpi[i][j] = camera.Xpi[i][j];
		}
	}
	return GZ_SUCCESS;
}


int GzRender::GzPushMatrix(GzMatrix	matrix)
{
	/* HW 3.9
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/
	/* checking for stack overflow */
	if (matlevel >= MATLEVELS) {
		return GZ_FAILURE;
	}
	/* matrix multiplication (dot-product) */
	if (matlevel <= 1) {
		/* identity matrix */
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				if (i == j) {
					Xnorm[normlevel][i][j] = 1;
				}
				else {
					Xnorm[normlevel][i][j] = 0;
				}
			}
		}
	}
	else {
		GzMatrix pureMatrix;
		makePureRotationMatrix(matrix, pureMatrix);
		dotProduct(pureMatrix, Xnorm, normlevel);

	}
	dotProduct(matrix, Ximage, matlevel);
	normlevel += 1;
	matlevel += 1;

	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
	/* HW 3.10
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/
	/* checking for stack overflow */
	if (matlevel <= 0) {
		return GZ_FAILURE;
	}
	matlevel -= 1;

	normlevel -= 1;
	return GZ_SUCCESS;
}





