#pragma once

#include<SFML/Graphics.hpp>

#include <fstream>
#include <filesystem>
//#include<string.h>
#include <sstream>


//#include<iostream>

#include"nnN1.h"
#include"Rocket.h"
#include"Graph.h"

#define rocCount 1300
#define drawRocCount 100
#define aMax 0.5
#define tresholdScore 1*pow(10,14)
#define reinitializetion 1500
#define itCount 100
#define lowTreshold 1*pow(10,7)
#define randCoef 0.01f

#define LoadNN 148
#define Fprogramm	2	//1 - train  2 - loadNN for play	3 - loadNN for train

#define DoMutate 1//1 - Y  0 - N
#define MakeBBS  1//1 - Y  0 - N

#define newParOnITBool		false
#define newParOnBESTBool	false
#define reinitBool			true


extern Graph graph;

using namespace std;
using std::filesystem::directory_iterator;
namespace fs = std::filesystem;

class train
{
public:








	Rocket  rocket[rocCount];


	sf::Vector3f roa;
	sf::Vector3f v;
	sf::Vector3f w;

	sf::Vector3f rOA[rocCount];
	sf::Vector3f rAB[rocCount];
	sf::Vector3f parentCol;

	sf::Color grColor;

	double aS;
	double bS;
	double gS;

	long double
		scorePar1 = 0,
		scorePar2 = 0,
		scorePar=0;

	sf::Vector3f col[rocCount];

	nnN1 baltazar[rocCount];
	nnN1 newnnN1[rocCount];
	nnN1 nnN1Best;
	nnN1 parent1, parent2;
	nnN1 baby;

	int gen = 0,
		it = 0,
		reinitIt=0,
		grLineN;

	ofstream outFile, outFileNN, outFileNNRead;
	ifstream inFileNN;

	float rocVisible[rocCount];
	int nBest;
	long double scoreSumm = 0;


	bool resKey = false;
	bool vis = true;













	void newParam() {
		//roa = sf::Vector3f{ float((rand() % 400) - 200),float((rand() % 400) - 200),-float((rand() % 200) - 50) } + sf::Vector3f(0.0f, 0.0f,-100.0f);
		//v = sf::Vector3f{ float(rand() % 50+10),float(rand() % 40 - 20),float(rand() % 40 - 20) };
		w = sf::Vector3f{ float(rand() % 100 / 100.0f * randCoef - 0.5f),float(rand() % 100 / 100.0f * randCoef - 0.5f),float(rand() % 100 / 100.0f * randCoef - 0.5f) };// float(rand() % 100 / 100.0f - 0.5f);;
		aS = double(rand() % 700 / 1000.0f * randCoef - 0.35f);
		bS = double(rand() % 700 / 1000.0f * randCoef - 0.35f);
		gS = double(rand() % 700 / 1000.0f * randCoef - 0.35f);
		//printf("\n\n\nx= %4f\ty= %4f\tvy= %4f\tw= %4f\ta= %4f \n ", roa.x, roa.y, v.y, w, a);

		roa = sf::Vector3f{abs( float((rand() % 400) - 200)),abs(float((rand() % 400) - 200)),-float((rand() % 200) - 50) } + sf::Vector3f(0.0f, 0.0f, -100.0f);
		v = sf::Vector3f{ float(rand() % 5 + 10),abs(float(rand() % 40 * randCoef - 20)),abs(float(rand() % 40 * randCoef - 20)) };





	}




	void initialize() {

		//testnnN1();

		newParam();
		printf("\n\n ============Params=============\nrOA {%2f ; %2f ; %2f}\tv {%2f ; %2f ; %2f}\n==============================\n", roa.x, roa.y, roa.z, v.x, v.y, v.z);
		//outFile.open("Scores.txt");

		setColor(rocCount, col);
		setRocket();
		if (Fprogramm == 1)
			for (int i = 0; i < rocCount; i++) {
				baltazar[i].initialize();
			}
		else if (Fprogramm == 2) {
			baltazar[0].initialize();
			fileInNN();
		}
		else if (Fprogramm == 3) {
			baltazar[0].initialize();
			fileInNN();
			for (int i = 1; i < rocCount; i++) {
				baltazar[i] = baltazar[0];
			}
		}


		
		
		//printf("\n\t\t________________%5f________________\n", baltazar[1].l[1].nodes[1].w[1]);
		//printf("\n\t\t%5f\n", adam[0].input.w[1][1]);
	}


	void setColor(int count, sf::Vector3f color[]) {
		color[0] = sf::Vector3f{ 1.0f,0.0f,0.0f };
		
		for (int i = 1; i < rocCount ; i++)color[i] = sf::Vector3f{ float(rand() % 100) / 100.0f,float(rand() % 100) / 100.0f ,float(rand() % 100) / 100.0f };
	}

	void setRocket() {


		for (int i = 0; i < rocCount; i++) {

			//rock[i].rOA = sf::Vector2f{ float((rand() % 100) - 50),float(rand() % 100 - 50) };
			rocket[i].rOA = roa;
			rocket[i].vm = v;
			rocket[i].wm = w;
			rocket[i].alpha = aS;
			rocket[i].betta = bS;
			rocket[i].gamma = gS;
			rocket[i].dead = false;
			rocket[i].cfuel = 0;
			rocket[i].tfuel = 0;
			rocket[i].engIgn = 1;
			rocket[i].engRun = false;
			rocket[i].initializeTurn();

		}

	}

	void play(sf::Shader& shader) {

		for (int i = 0; i < rocCount; i++) {
			if (!rocket[i].dead) {
				rocket[i].dynamic();
				//rOA[i] = rocket[i].rOA;
				//rAB[i] = rocket[i].rAB;
				//rocVisible[i] = rocket[i].visible;
			}

		}
		
		rocket[0].traektory(shader);
		drawAll(shader);
		//rocket[0].draw(shader);
		//printf("%10f\t%1f\t%1f\t%3f\n", rocket[0].rOA.x, rocket[0].rOA.y, rocket[0].rOA.z,1.0f );

		if (allDead()) {

			showBest();
			if(MakeBBS==1) makeBabies();
			
			if(DoMutate==1)mutateBabies();
			

			restart();

		}

	}


	void play_NoVisual() {
	
		for (int i = 0; i < rocCount; i++) {
			if (!rocket[i].dead) {
				rocket[i].dynamic();
			}
		}
		if (allDead()) {
			printf("\nNV");
			showBest();
			printf("\n");
			if (MakeBBS == 1) makeBabies();

			if (DoMutate == 1)mutateBabies();
			restart();
		}
	}

	bool allDead() {

	for (int i = 0; i < rocCount; i++) {

		if ((-rocket[i].rOA.z <= 0 | -rocket[i].rOA.z > 500 | abs(rocket[i].rOA.x) > 500 | rocket[i].rOA.y > 500) & !rocket[i].dead) {

			rocket[i].dead = true;
			baltazar[i].score = calcScore(rocket[i]);

		}
	}

for (int i = 0; i < rocCount; i++) {
	if (!rocket[i].dead)
		return(false);
}

return(true);

	}

	double calcScore(Rocket rock) {
		//if (y == 0)y = 0.01;
		float y = rock.rOA.y;
		float x = rock.rOA.x;
		float z = rock.rOA.z;
		float a1 = rock.alpha;
		float b1 = rock.betta;
		float g1 = rock.gamma;
		double cf = rock.cfuel;
		double tf = rock.tfuel;
		float vy = rock.vm.y;
		float vx = rock.vm.x;
		float vz = rock.vm.z;
		float rDot =100* rock.dot(rock.normalize(rock.vm), rock.normalize(rock.vecX));

		long double score =
			(pow(2, -pow(z, 2) / 100) * 5 + pow(2, -pow(z, 2) / 5000) * 1) *
				(
					(
						(
							(pow(2, -pow(vz, 2) / 250) * 50 + (pow(2, -pow(vz, 2) * 2) * 50)) *
							(pow(2, -pow(vx, 2) / 150) * 200 + (pow(2, -pow(vx, 2) * 4) * 200)) *
							(pow(2, -pow(vy, 2) / 250) * 50 + (pow(2, -pow(vy, 2) * 2) * 50)) 
						)*
						(atan(-pow(tf - 202, 2) + 200 * 200) + 1.57)+
						

						(
							(pow(2, -pow(b1, 2) * 10) * 50 + pow(2, -pow(b1, 2) / 10) * 10) *
							(pow(2, -pow(g1, 2) * 10) * 50 + pow(2, -pow(g1, 2) / 10) * 10) *
							(pow(2, -pow(a1, 2) * 4) * 25 + pow(2, -pow(a1, 2) / 10) * 10) 
							
						)*
						
						(atan(-pow(cf - 200, 2) + 200 * 200) + 1.6) *
						(pow(2, -pow(rDot+100, 2) /70) * 2000 + pow(2, -pow(rDot+100, 2) / 1000) * 5)
						
					+

					(
						(pow(2, -pow(vz, 2) / 200) * 20 + (pow(2, -pow(vz, 2) * 2) * 300)) *
						(pow(2, -pow(vx, 2) / 50) * 100 + (pow(2, -pow(vx, 2) *4) * 15000)) *
						(pow(2, -pow(vy, 2) / 200) * 20 + (pow(2, -pow(vy, 2) * 2) * 300)) *


						(pow(2, -pow(b1, 2) * 10) * 50 + pow(2, -pow(b1, 2) / 10) * 2) *
						(pow(2, -pow(g1, 2) * 10) * 50 + pow(2, -pow(g1, 2) / 10) * 2) *
						(pow(2, -pow(a1, 2) * 4) * 25 + pow(2, -pow(a1, 2) / 10) * 2) *

						(atan(-pow(cf - 200, 2) + 200 * 200) +1.572)*
						(atan(-pow(tf - 202, 2) + 200 * 200) + 1.57)*
						(pow(2, -pow(rDot + 100, 2) / 70) * 2000 + pow(2, -pow(rDot + 100, 2) / 1000) * 10)


					)


					)

				) / ( rocCount) + 0.001;		//* ((pow(2, -pow(cf, 2) /5000) * 50 + (pow(2, -pow(cf, 2) /50) * 200))				//(pow(2, -pow(y, 2) / 50) * 400 + 50 - abs(0.1f *y)) * (pow(2, -pow(a, 2) /2)*100+50 -abs(0.25f*a))  * (pow(2, -pow(vy, 2) / 2) * 80 + 10 - abs(0.05f * vy))						(pow(2, -pow(y, 2) / 50) * 400 + 50 - abs(0.1f *y))* (pow(2, -pow(vy, 2) / 2) * 100 + 50 - abs(0.25f *vy))
				
		if (score < 0)score = 0.001;

		







		return(score);
	}

	void makeBabies() {
		calcSumScore();
		printf("\n___%f___\n", scoreSumm);
		findBest();
		
	
		baltazar[0] = baltazar[nBest];
		col[0] = col[nBest];
		int parint;


		for (int i = 1; i < rocCount; i++) {
			//newnnN1[i] = baltazar[selectrParent()];
			newnnN1[i] = sexBreeding(baltazar[selectrParent()], baltazar[selectrParent()]);
		}
		for (int i = 1; i < rocCount; i++) {
			baltazar[i] = newnnN1[i];
		}

		//sf::Vector3f col1, col2;
		/*
		for (int i = 1; i < rocCount; i++) {

			parint = selectrParent();
			
			for (int lnum = 1; lnum <= layerCount; lnum++) {
				for (int j = 1; j <= parent1.l[lnum].nCount; j++) {
					parent1.l[lnum].nodes[j] = baltazar[0].l[lnum].nodes[j];
				}
			}
			
			//nnCopy(baltazar[parint], parent1);
			scorePar1 = scorePar;
			col1 = parentCol;
			
			parint = selectrParent();

			
			for (int lnum = 1; lnum <= layerCount; lnum++) {
				for (int j = 1; j <= parent1.l[lnum].nCount; j++) {
					parent2.l[lnum].nodes[j] = baltazar[0].l[lnum].nodes[j];
				}
			}
			
			//nnCopy(baltazar[parint], parent2);
			scorePar2 = scorePar;
			col2 = parentCol;
			
			if (rand() % 2 == 1) col[i] = col1;
			else  col[i] = col2;
			
			//nnSexCopy(baltazar[selectrParent()], baltazar[selectrParent()], newnnN1[i]);

		
		
		
		for (int lnum = 2; lnum <= layerCount; lnum++) {
				for (int j = 1; j <= newnnN1[i].l[lnum].nCount; j++) {
					if (rand() % 2 == 1) newnnN1[i].l[lnum].nodes[j] = parent1.l[lnum].nodes[j];
					else newnnN1[i].l[lnum].nodes[j] = parent2.l[lnum].nodes[j];
				}
			}
		}
		for (int i = 1; i < rocCount; i++) {
			for (int lnum = 2; lnum <= layerCount; lnum++) {
				for (int j = 1; j <= newnnN1[i].l[lnum].nCount; j++) {
					 baltazar[i].l[lnum].nodes[j] = newnnN1[i].l[lnum].nodes[j];
					
				}
			}
			
		}
		*/
		
		
				

		
		
		
	}

	nnN1 sexBreeding(nnN1 p1, nnN1 p2) {
		nnN1 fetus = p1;
		for (int lnum = 2; lnum <= layerCount; lnum++) {
			for (int j = 1; j <= fetus.l[lnum].nCount; j++) {
				if (rand() % 2 == 1) fetus.l[lnum].nodes[j] = p2.l[lnum].nodes[j];
			}
		}
		return(fetus);
	
	}

	void calcSumScore() {
		scoreSumm = 0;

		for (int i = 0; i < rocCount; i++)
			scoreSumm += baltazar[i].score;


	}

	void findBest() {
		long double bestScore = baltazar[0].score;
		nBest = 0;
		for (int i = 1; i < rocCount; i++) {
			if (bestScore < baltazar[i].score) {
				bestScore = baltazar[i].score;
				nBest = i;
			}
		}

	}

	int selectrParent() {
		long double runSum = 0;
		
		long double randomScore =  rand() % (long int)round(scoreSumm);
		
		for (int i = 0; i < rocCount  ; i++) {
			runSum += baltazar[i].score;
			if (runSum >= randomScore ) {
				parentCol = col[i];
				scorePar = baltazar[i].score;
				
				return(i); 
					
			}
		}
	}

	

	

	void mutateBabies() {

		for (int i = 1; i < rocCount; i++) {

			baltazar[i].mutation();


			float red = float(rand() % 20) / 100.0f - 0.1f;
			float green = float(rand() % 20) / 100.0f - 0.1f;
			float blue = float(rand() % 20) / 100.0f - 0.1f;
			if ((col[i].x + red) < 0.9f & (col[i].x + red) > 0.2f)col[i].x += red;
			if ((col[i].y + green) < 0.9f & (col[i].y + green) > 0.2f)col[i].y += green;
			if ((col[i].z + blue) < 0.9f & (col[i].z + blue) > 0.2f)col[i].z += blue;

			if ((col[i].x * col[i].x + col[i].y * col[i].y + col[i].z * col[i].z < 0.008f) | (col[i].x * col[i].x + col[i].y * col[i].y + col[i].z * col[i].z > 0.73f))col[i] = sf::Vector3f{ float(rand() % 100) / 100.0f,float(rand() % 100) / 100.0f ,float(rand() % 100) / 100.0f };

		}
	}


	void restart() {

		resKey = false;
		
		//printf("\n\t\t________________%5f________________\n", adam[0].input.w[1][1]);

		//addScore(nnN1Best.score);

		if (bdrawGraph) {
			graph.drawline(baltazar[0].score, grLineN,grColor);
			graph.lines[grLineN].x += 1;
			//graph.lines[0].text.setString(std::to_string(eMax));
			graph.lines[grLineN].bdraw = true;
		}



		if (it > itCount) {
			if (baltazar[0].score > tresholdScore) {
				system("cls");


				gen++;
				showBest();
				//printf("_________________________________Generation %d______________________________", gen);

				

				fileOutNN();
				fileOutNN_forRead();
				if (newParOnBESTBool) {
					newParam();
					//it = 0;
				}

			}
			if (newParOnITBool) {
				newParam();
				//it = 0;
			}
			it = 0;
		}
		it++;
		
		if (reinitIt >= reinitializetion & reinitBool) {
			if (baltazar[0].score > lowTreshold) reinitIt = 0;
			else { 
				reinitIt = 0;
				gen = 0;
				for (int i = 0; i < rocCount; i++) {
					baltazar[i].initialize();
				}
			}
		}
		
		reinitIt++;
		setRocket();
	}

	void restart_NoVisual() {

		


		resKey = false;
		
		if (it > itCount) {
			if (baltazar[0].score > tresholdScore) {
				

				fileOutNN();
				fileOutNN_forRead();
				if (newParOnBESTBool) {
					newParam();
					it = 0;
				}
			}
			if (newParOnITBool) {
				newParam();
				it = 0;
			}
		}
		it++;
		
		if (reinitIt >= reinitializetion & reinitBool) {
			if (baltazar[0].score > lowTreshold) reinitIt = 0;
			else {
				reinitIt = 0;
				gen = 0;
				for (int i = 0; i < rocCount; i++) {
					baltazar[i].initialize();
				}
			}
		}
		
		reinitIt++;

		setRocket();
		
	}

	void showBest() {

		double bestScore = baltazar[0].score;
		float xbest = rocket[0].rOA.x;
		float ybest = rocket[0].rOA.y;
		float zbest = rocket[0].rOA.z;
		double abest = rocket[0].alpha;
		double bbest = rocket[0].betta;
		double gbest = rocket[0].gamma;
		float vybest = rocket[0].vm.y;
		float vxbest = rocket[0].vm.x;
		float vzbest = rocket[0].vm.z;
		int engIgnB = rocket[0].engIgn;
		double cf = rocket[0].cfuel;
		double tf = rocket[0].tfuel;
		float rDot = rocket[0].dot(rocket[0].normalize(rocket[0].vm), rocket[0].normalize(rocket[0].vecX));
		/*
		for (int i = 1; i < rocCount; i++) {
			if (bestScore < baltazar[i].score) {
				 xbest = rocket[i].rOA.x;
				 ybest = rocket[i].rOA.y;
				 zbest = rocket[i].rOA.z;
				 abest = rocket[i].alpha;
				 bbest = rocket[i].betta;
				 gbest = rocket[i].gamma;
				 vybest = rocket[i].vm.y;
				 vxbest = rocket[i].vm.x;
				 vzbest = rocket[i].vm.z;
				 engIgnB = rocket[i].engIgn;
				 cf = rocket[i].cfuel;
				

			}
		}
		*/
		printf("\n reInit %d", it);
		printf("\nbest score %5f  a=%4f  b=%4f  g=%4f  vx=%4f vy=%4f  vz=%4f  cf=%2f  tf=%2f  ign=%d  dot=%f", bestScore, abest,bbest,gbest,vxbest,vybest,vzbest,cf,tf, engIgnB,rDot);
	}


	void control() {

		float inp[10], outp[4];

		for (int i = 0; i < rocCount; i++) {
			if (!rocket[i].dead) {
				inp[0] = rocket[i].vm.x;
				inp[1] =( rocket[i].vm.y);
				inp[2] =( rocket[i].vm.z);
				inp[3] = rocket[i].wm.x;
				inp[4] = rocket[i].wm.y;
				inp[5] = rocket[i].wm.z;
				inp[6] = rocket[i].alpha;
				inp[7] = rocket[i].betta;
				inp[8] = rocket[i].gamma;
				inp[9] = -rocket[i].rOA.z;
				
				baltazar[i].thinK(inp, outp);

				//outp[0] - thrust	bthrust
				//outp[1] - roll	bturnX
				//outp[2] - pitch	bturnY
				//outp[3] - yaw		bturnZ

				printf("\n == %f \t %f \t %f\t %f", outp[0], outp[1], outp[2], outp[3]);
				// bthrust
				if (outp[0] < 0.9f)outp[0] = 0;
				else {
					if (rocket[i].engRun)
						rocket[i].tfuel+=1* outp[0];
				}
				if (rocket[i].tfuel > 600)outp[0] = 0;


				// bturnX
				if (outp[1] > 0.7f) {
					outp[1] = 1;
					
					rocket[i].cfuel += 1;
				}
				else if (outp[1] < -0.7f) {		/////////////////////////////////////////////////////////// sigmoid -1 .. +1
					outp[1] = -1;
					rocket[i].cfuel += 1;
				}
				else outp[1] = 0;

				// bturnY
				if (outp[2] > 0.6f) {
					outp[2] = 1;
					rocket[i].cfuel += 1;
				}
				else if (outp[2] < -0.6f) {		/////////////////////////////////////////////////////////// sigmoid -1 .. +1
					outp[2] = -1;
					rocket[i].cfuel += 1;
				}
				else outp[2] = 0;

				// bturnZ
				if (outp[3] > 0.6f) {
					outp[3] = 1;
					rocket[i].cfuel += 1;
				}
				else if (outp[3] < -0.6f) {		/////////////////////////////////////////////////////////// sigmoid -1 .. +1
					outp[3] = -1;
					rocket[i].cfuel += 1;
				}
				else outp[3] = 0;

				rocket[i].bThrust = outp[0];
				rocket[i].bTurnX = outp[1];
				rocket[i].bTurnY = outp[2];
				rocket[i].bTurnZ = outp[3];
				
				rocket[i].control();

			}
		}

	}

	void drawAll(sf::Shader& shader) {

	
		sf::Vector3f urOA[500];
		sf::Vector3f urOB[500];
		float urockCount = drawRocCount;
		sf::Vector3f uCol[500];
		
		
	





		for (int i = 0; i < drawRocCount; i++) {
			rocket[i].preDraw();

			urOA[i] = rocket[i].rOA;
			urOB[i] = rocket[i].rOB;
			uCol[i] = col[i];



		}

		shader.setUniformArray("u_coneROA", urOA,urockCount);
		shader.setUniformArray("u_coneROB", urOB,urockCount);
		shader.setUniformArray("u_rocCol", uCol, urockCount);
		shader.setUniform("u_rockCount", urockCount);

		shader.setUniform("u_thrust", rocket[0].vecThrust);
		shader.setUniform("u_lvecCtlr", rocket[0].lvecCtlr);
		shader.setUniformArray("u_vecCtrl", rocket[0].vecCtrl, 7);
		shader.setUniformArray("u_fire", rocket[0].fireDraw, 7);
		shader.setUniform("u_ra", rocket[0].ra);
		shader.setUniform("u_rb", rocket[0].rb);


	}


	void addScore(double score1) {

		/*

		//outFile.open("weights.txt", ios::app);

		if (outFile.is_open()) {
			outFile << score1 << "\n";

		}
		else {
			outFile.open("Scores.txt");
			outFile << score1 << "\n";
		}

		*/


	}
	
	void fileOutNN_forRead() {


		string path = dirPath() + "\\Baltazar\\NNFiles_Read";

		int i = 0;
		
		directory_iterator end_itr;
		for (directory_iterator itr(path); itr != end_itr; ++itr)
		{
			fs::path path1{ itr->path() };
			std::string path_string{ path1.string() };
			i++;

		}
		
		std::stringstream ss;
		ss << path << "\\NNFileRead_" << i << ".txt";
		std::string str = ss.str();
		ss.str("");


	
		
		outFileNNRead.open(str);
		outFileNNRead
			<< "\n======================================================"
			<< "\nScore\t" << baltazar[0].score
			<< "\nZ\t" << rocket[0].rOA.z
			<< "\nalpha\t" << rocket[0].alpha
			<< "\nbetta\t" << rocket[0].betta
			<< "\ngamma\t" << rocket[0].gamma
			<< "\nVx\t" << rocket[0].vm.x
			<< "\nVy\t" << rocket[0].vm.y
			<< "\nVz\t" << rocket[0].vm.z
			<< "\nengIng\t" << rocket[0].engIgn
			<< "\ncFuel\t" << rocket[0].cfuel
			<< "\ntfuel\t" << rocket[0].tfuel
			<< "\n======================================================"
			<< "\n\nInitial Parametrs"
			<< "\nrOa\t" << "{ " << roa.x << " ; " << roa.y << " ; " << roa.z << " }\n"
			<< "\nV\t" << "{ " << v.x << " ; " << v.y << " ; " << v.z << " }\n"
			<< "\nw\t" << "{ " << w.x << " ; " << w.y << " ; " << w.z << " }\n"
			<< "\nalpha\t" << aS
			<< "\nbetta\t" << bS
			<< "\ngamma\t" << gS
			<< "\n======================================================";
		int nEnableK = 0;
		/*
		for (int lnum = 1; lnum <= layerCount; lnum++) {
			outFileNNRead << "\n\n Layer " << lnum << "\n";
			for (int i = 1; i <= baltazar[0].l[lnum].nCount; i++) {
				outFileNNRead << "\nnode  " << i << "\t\t link" << "\t weight" << "\t bias " << baltazar[0].l[lnum].nodes[i].bias << "\n";
				for (int k = 1; k <= nodeCount; k++) {
					nEnableK = 0;
					if (baltazar[0].l[lnum].nodes[i].n[k])nEnableK = 1;
					outFileNNRead << "\t\t" << k << "  " << nEnableK << "\t" << baltazar[0].l[lnum].nodes[i].w[k] << "\n";
				}
			}
		}
		*/
		for (int lnum = 2; lnum <= layerCount; lnum++) {
			outFileNNRead << "\n\n Layer " << lnum << "\n";
			for (int i = 1; i <=baltazar[0].l[lnum].nCount; i++) {
				outFileNNRead << "\nnode  " << i << "\t\t link" << "\t weight" << "\t bias " << baltazar[0].l[lnum].nodes[i].bias << "\n";
				for (int j = 1; j <= (lnum - 1); j++) {
					for (int k = 1; k <= baltazar[0].l[j].nCount; k++) {
						nEnableK = 0;
						if (baltazar[0].l[lnum].nodes[i].n[baltazar[0].l[j].nodes[k].nNum])nEnableK = 1;
						outFileNNRead << "\t\t" << baltazar[0].l[j].nodes[k].nNum << "  " << nEnableK << "\t" << baltazar[0].l[lnum].nodes[i].w[baltazar[0].l[j].nodes[k].nNum] << "\n";
							
					}
				}
				
			}
		}
		
		outFileNNRead.close();
		
	}

	void fileOutNN() {


		string path = dirPath() + "\\Baltazar\\NNFiles";
		
		int i = 0;
		
		directory_iterator end_itr;
		for (directory_iterator itr(path); itr != end_itr; ++itr)
		{
			fs::path path1{ itr->path() };
			std::string path_string{ path1.string() };
			i++;

		}
		
		std::stringstream ss;
		ss << path<<"\\NNFile_" << i << ".txt";
		std::string str = ss.str();
		ss.str("");


		outFileNN.open(str);

		
		int nEnableK = 0;
		for (int lnum = 1; lnum <= layerCount; lnum++) {
			outFileNN << "\n\n" ;
			for (int i = 1; i <= baltazar[0].l[lnum].nCount; i++) {
				outFileNN << baltazar[0].l[lnum].nodes[i].bias << "\n";
				for (int k = 1; k <= nodeCount; k++) {
					nEnableK = 0;
					if (baltazar[0].l[lnum].nodes[i].n[k])nEnableK = 1;
					outFileNN << nEnableK << "\n";
					outFileNN << baltazar[0].l[lnum].nodes[i].w[k] << "\n";
				}
			}
		}

		outFileNN.close();

	}

	void fileInNN() {

		string path = dirPath() + "\\Baltazar\\NNFiles";
		std::stringstream ss;
		ss << path << "\\NNFile_" << LoadNN << ".txt";
		std::string str = ss.str();
		ss.str("");

		inFileNN.open(str);

		int nEnableK = 0;

		for (int lnum = 1; lnum <= layerCount; lnum++) {
			
			for (int i = 1; i <= baltazar[0].l[lnum].nCount; i++) {
				inFileNN>>baltazar[0].l[lnum].nodes[i].bias  ;
				for (int k = 1; k <= nodeCount; k++) {
					inFileNN >> nEnableK;
					if (nEnableK == 1)baltazar[0].l[lnum].nodes[i].n[k] = true;
					else baltazar[0].l[lnum].nodes[i].n[k] = false;
					inFileNN >> baltazar[0].l[lnum].nodes[i].w[k];
				}
			}
		}

		inFileNN.close();

	}

	void fileOutRez(int fnum) {
		string path = dirPath() + "\\Baltazar\\Reserve";

		

		std::stringstream ss;
		ss << path << "\\NNFile_Rez" << fnum << ".txt";
		std::string str = ss.str();
		ss.str("");


		outFileNN.open(str);


		int nEnableK = 0;
		for (int lnum = 1; lnum <= layerCount; lnum++) {
			outFileNN << "\n\n";
			for (int i = 1; i <= baltazar[0].l[lnum].nCount; i++) {
				outFileNN << baltazar[0].l[lnum].nodes[i].bias << "\n";
				for (int k = 1; k <= nodeCount; k++) {
					nEnableK = 0;
					if (baltazar[0].l[lnum].nodes[i].n[k])nEnableK = 1;
					outFileNN << nEnableK << "\n";
					outFileNN << baltazar[0].l[lnum].nodes[i].w[k] << "\n";
				}
			}
		}

		outFileNN.close();
	}


	string dirPath() {
		/*
		TCHAR buffer[MAX_PATH] = { 0 };
		GetModuleFileName(NULL, buffer, MAX_PATH);
		std::wstring::size_type pos = std::wstring(buffer).find(L'x64');
		//wstring pos1 = std::wstring(buffer).append(L"\\florist\\trainSamples\\");
		wstring wstr = std::wstring(buffer).substr(0, pos);
		string str(wstr.begin(), wstr.end());;
		*/
		return ("C:\\Users\\themr\\Desktop\\New folder 2\\VS\\Baltazar");
	}
	
};

