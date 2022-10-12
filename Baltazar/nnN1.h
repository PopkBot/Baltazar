#pragma once
#include<thread>
#include <stdio.h>    
#include <stdlib.h> 

#include<iostream>
#include <chrono>
#include<time.h>

#include<math.h>
#include"Rocket.h"


#define xfac 10.0f
#define mutationLinkRate 0.02	//<1
#define mutationWandBRate 0.1	//<1


#define inputNodesCount 10
#define outputNodesCount 4
//#define h1NodesCount 6
//#define h2NodesCount 5

#define hiddenLayersCount 3
#define hiddenNodesCount 10
#define nodeCount inputNodesCount+outputNodesCount+hiddenLayersCount*hiddenNodesCount+1

#define layerCount 2+hiddenLayersCount

#define mutateLinks true
#define mutateWandB true

class nnN1
{
public:
	struct node
	{
		bool n[nodeCount];
		double w[nodeCount];
		double bias;
		double inVal;
		double outVal;
		bool nEnable;
		int nNum;
	};

	struct layer
	{
		node nodes[15];
		int nCount;

	};

	float linksCount = 0;
	float enableNodeCount;
	float maxLinks;
	float maxEnNodes = nodeCount;

	long double score;



	layer l[layerCount + 1];

	node nodes[nodeCount + layerCount];

	void initialize() {
		layer input,output, h[hiddenLayersCount];
		//srand(time(0));

		for (int i = 0; i <= nodeCount; i++) {
			nodes[i].nNum = i;
			for (int j = 0; j <= nodeCount; j++)
			{
				nodes[i].n[j] = false;
				nodes[i].w[j] = 0;// rand() % 2000 / 100.0f - 10.0f;
				nodes[i].bias = 0;
				nodes[i].inVal = 0;
				nodes[i].outVal = 0;
			}
		}

		int nodeN = 1;

		input.nCount = inputNodesCount;
		output.nCount = outputNodesCount;

		for (int i = 1; i <= inputNodesCount; i++) {
			input.nodes[i] = nodes[nodeN];
			nodeN++;
		}
		for (int i = 1; i <= outputNodesCount; i++) {
			output.nodes[i] = nodes[nodeN];
			nodeN++;
		}
		for (int j = 0; j < hiddenLayersCount; j++) {
			for (int i = 1; i <= hiddenNodesCount; i++) {
				h[j].nCount = hiddenNodesCount;
				h[j].nodes[i] = nodes[nodeN];
				nodeN++;
			}
		}


		l[1] = input;
		//l[2] = h[0];
		//l[3] = h[1];
		l[layerCount] = output;

		for (int i = 2; i < layerCount; i++) {
			l[i] = h[i - 2];
		}


		
		/*
		int nodeN = 1;
		for (int i = 1; i <= inputNodesCount; i++) {
			input.nodes[i] = nodes[nodeN];
			nodeN++;
		}
		for (int i = 1; i <= outputNodesCount; i++) {
			output.nodes[i] = nodes[nodeN];
			nodeN++;
		}
		for (int i = 1; i <= h1NodesCount; i++) {
			h1.nodes[i] = nodes[nodeN];
			nodeN++;
		}
		for (int i = 1; i <= h2NodesCount; i++) {
			h2.nodes[i] = nodes[nodeN];
			nodeN++;
		}

		l[1] = input;
		l[2] = h1;
		l[3] = h2;
		l[4] = output;
		*/




		

		for (int i = 2; i <= layerCount; i++) {
			maxLinks += l[i - 1].nCount * l[i].nCount;
		}

		printf("\n\n= = = = = = NN initialized = = = = = =\n\n");
	}




	float sigmoid(float val) {
		
		return(1 / (1 + 1 / exp(val)) * 2 - 1);
		
	}

	float x_10(float val) {
		float func;

		func = val / 10;
		if (func > 1)func = 1;
		if (func < -1)func = -1;
		
		return(func);
	}

	void thinK(float inp[], float outp[]) {

		for (int i = 1; i <= l[1].nCount; i++) {
			l[1].nodes[i].outVal = (double)inp[i - 1];
		}

		layersCalc();

		for (int i = 1; i <= l[layerCount].nCount; i++) {

			outp[i - 1] = (float)l[layerCount].nodes[i].outVal;
		}
		//printf("\n\n=== %f\n\n", l[layerCount].nodes[outputNodesCount].outVal);
	}

	void layersCalc() {
		double val = 0;
		for (int lnum = 2; lnum <= layerCount; lnum++) {
			for (int i = 1; i <= l[lnum].nCount; i++) {
				val = 0;
				for (int j = 1; j <= (lnum - 1); j++) {
					for (int k = 1; k <= l[j].nCount; k++) {
						if (l[lnum].nodes[i].n[l[j].nodes[k].nNum])
							val += l[j].nodes[k].outVal * l[lnum].nodes[i].w[l[j].nodes[k].nNum];
					}
				}
				val += l[lnum].nodes[i].bias;
				l[lnum].nodes[i].outVal = sigmoid(val);
			}
		}


	}

	void mutation() {

		int bn;
		//if ((rand() % 20000 / 1000.0f < mutationLinkRate) & mutateLinks) initialize();

		for (int lnum = 2; lnum <= layerCount; lnum++) {
			for (int i = 1; i <= l[lnum].nCount; i++) {
				if ((rand() % 10000 / 10000.0f < mutationLinkRate) & mutateLinks) {
					bn = rand() % (nodeCount + 1);
					l[lnum].nodes[i].n[bn] = !l[lnum].nodes[i].n[bn];
				}
				if ((rand() % 1000 / 1000.0f < mutationWandBRate) & mutateWandB) {
					bn = rand() % (nodeCount + 1);
					l[lnum].nodes[i].w[bn] += rand() % (int)(xfac * 1000) / 1000.0f - xfac / 2.0f;
				}
				if ((rand() % 1000 / 1000.0f < mutationWandBRate) & mutateWandB)
					l[lnum].nodes[i].bias += rand() % (int)(xfac * 1000) / 1000.0f - xfac / 2.0f;


			}

			linksNodesCount();
		}



	}

	void linksNodesCount() {

		for (int i = 1; i <= l[1].nCount; i++) l[1].nodes[i].nEnable = true;
		for (int i = 1; i <= l[layerCount].nCount; i++) l[layerCount].nodes[i].nEnable = true;

		for (int lnum = 2; lnum < layerCount; lnum++) {
			for (int i = 1; i <= l[2].nCount; i++) l[lnum].nodes[i].nEnable = false;
		}
		
		enableNodeCount = l[1].nCount + l[layerCount].nCount;
		linksCount = 0;

		for (int lnum = 2; lnum <= layerCount; lnum++) {
			for (int i = 1; i <= l[lnum].nCount; i++) {
				for (int j = 1; j <= (lnum - 1); j++) {
					for (int k = 1; k <= l[j].nCount; k++) {
						if (l[lnum].nodes[i].n[l[j].nodes[k].nNum]) {
							l[lnum].nodes[i].nEnable = true;
							l[j].nodes[k].nEnable = true;
							linksCount++;
						}
					}
				}
			}
		}
		for (int lnum = 2; lnum <= layerCount - 1; lnum++) {
			for (int i = 1; i <= l[lnum].nCount; i++) {
				if (l[lnum].nodes[i].nEnable)
					enableNodeCount++;
			}
		}

	}



	float sigmDeriv(float val) {

		return(val * (1 - val));
	}

};

