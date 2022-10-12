#pragma once


#include <SFML/Graphics.hpp>
#include<iostream>
#include <algorithm>
#include <vector>
#include <math.h>

#define wg 1000
#define hg 500
#define lineCount 4

#define maxH 10
//#define TRHly hg/3

#define bdrawGraph true
#define gridNum 21

using namespace std;
class Graph
{

public:

	sf::Shader shader;

	struct StLines {
		std::vector<sf::Vertex>vertex;
		float x;
		bool bDrawLine = false;
		sf::Text text;

		sf::Vertex posLine_Y[2];
		bool bdraw = false;
	};

	
	struct StMarkers {
		sf::CircleShape cShape = sf::CircleShape(10);
		sf::Text markerText;
		bool mdraw = false;


	};

	StMarkers marker[lineCount];

	sf::Font font;

	struct Axis {
		std::vector<sf::Vertex>xAxis, yAxis;
		float x0;
		float y0;
	};

	Axis axis;
	Axis axisR;
	StLines grid[gridNum];

	sf::Text tmaxV;

	StLines lines[lineCount];

	StLines treashHoldLine;

	bool bCheck_y = false,
		bScale_y = false;
	float
		qy = 2,
		dqy = 0;


	void initialize() {
		sf::RenderWindow windowg(sf::VideoMode(wg, hg), "Tracing", sf::Style::Titlebar | sf::Style::Close); //Titlebar | sf::Style::Close
		windowg.setFramerateLimit(60);
		windowg.setPosition(sf::Vector2i(830, 0));

		sf::RenderTexture empityTexture;
		empityTexture.create(wg, hg);
		sf::Sprite empitySprite = sf::Sprite(empityTexture.getTexture());



		shader.setUniform("u_resolution", sf::Vector2f(wg, hg));


		if (!font.loadFromFile("arial.ttf")) {
			printf("\nerror=======================================\n");
		}


		tmaxV.setPosition(axisR.x0 + 10, 30);
		tmaxV.setFillColor(sf::Color{ 150,150,150 });
		tmaxV.setFont(font);
		tmaxV.setCharacterSize(15);



		sf::RenderStates states;
		sf::Transform transform;
		states.transform = transform;

		while (windowg.isOpen()) {


			sf::Event event;

			while (windowg.pollEvent(event))
			{



				if (event.type == sf::Event::Closed)
				{

					windowg.close();
				}

				if (event.type == sf::Event::KeyPressed) {
					if (event.key.code == sf::Keyboard::Escape) {
						windowg.close();
					}
					if (event.key.code == sf::Keyboard::Space) {
						//drawline(10, 0);
					}
					if (event.key.code == sf::Keyboard::W) {



					}
					if (event.key.code == sf::Keyboard::S) {




					}

				}
				else if (event.type == sf::Event::KeyReleased)
				{
					if (event.key.code == sf::Keyboard::Space) {

					}
					if (event.key.code == sf::Keyboard::W) {

					}


				}

			}


		//	tmaxV.setString(to_string(axis.y0 / (qy)));



			windowg.clear();
			//boundscheck();

			//check_y();
			//scale_y_max();

			windowg.draw(&axis.xAxis[0], axis.xAxis.size(), sf::LinesStrip);
			windowg.draw(&axis.yAxis[0], axis.yAxis.size(), sf::LinesStrip);
			windowg.draw(&axisR.yAxis[0], axisR.yAxis.size(), sf::LinesStrip);
		//	windowg.draw(tmaxV);
			windowg.draw(&treashHoldLine.posLine_Y[0], 2, sf::LinesStrip);
			windowg.draw(treashHoldLine.text);

			for (int i = 0; i < lineCount; i++) {
				if (lines[i].bdraw) {
				//	windowg.draw(lines[i].text);
				//	windowg.draw(&lines[i].posLine_Y[0], 2, sf::LinesStrip);
					windowg.draw(&lines[i].vertex[0], lines[i].vertex.size(), sf::LinesStrip);

					
				}
				if (marker[i].mdraw) {
					windowg.draw(marker[i].markerText);
					windowg.draw(marker[i].cShape);
				}

			}
			
			for (int i = 0; i < gridNum; i++) {
				windowg.draw(&grid[i].posLine_Y[0], 2, sf::LinesStrip);
				windowg.draw(grid[i].text);
			}


			windowg.display();










		}
	}

	void gridInitialize() {

		float step = axis.y0/ (gridNum - 1);

		qy = step;

		for (int i = 0; i < gridNum; i++) {
			grid[i].text.setString(to_string(i));
			grid[i].text.setPosition(axisR.x0 + 10, axis.y0-step*i);
			grid[i].text.setFillColor(sf::Color{ 150,150,150 });
			grid[i].text.setFont(font);
			grid[i].text.setCharacterSize(10);

			grid[i].posLine_Y[1].position.x = axisR.x0 + 5;
			grid[i].posLine_Y[1].position.y = axis.y0 - step*i;
			grid[i].posLine_Y[1].color = sf::Color{ 100,100,100 };

			grid[i].posLine_Y[0].position.x = axisR.x0 - 5;
			grid[i].posLine_Y[0].position.y = axis.y0 - step * i;
			grid[i].posLine_Y[0].color = sf::Color{ 100,100,100 };

		}

	}

	void markerInitialize(int mnum, sf::Color mcol,string mtext) {

		marker[mnum].cShape.setFillColor(mcol);
		marker[mnum].cShape.setPosition(axisR.x0+30,hg/2+25*mnum);
		marker[mnum].markerText.setString(mtext);
		marker[mnum].markerText.setPosition(axisR.x0 + 60, hg/2 + 25 * mnum);
		marker[mnum].markerText.setFillColor(sf::Color{ 150,150,150 });
		marker[mnum].markerText.setFont(font);
		marker[mnum].markerText.setCharacterSize(15);
		marker[mnum].mdraw = true;


	}

	void initializeTHLine(double trh) {
		

		//qy = (axis.y0 - TRHly) / log10f(trh);
		printf("\n\n=========qy = %f\n", qy);

		treashHoldLine.posLine_Y[1].position.y = axis.y0-qy*log10f(trh);
		treashHoldLine.posLine_Y[1].position.x = axisR.x0;
		treashHoldLine.posLine_Y[1].color = sf::Color{ 50,100,50 };

		treashHoldLine.posLine_Y[0].position.x = axis.x0;
		treashHoldLine.posLine_Y[0].position.y = axis.y0 - qy * log10f(trh);
		treashHoldLine.posLine_Y[0].color = sf::Color{ 50,100,50 };



		treashHoldLine.text.setPosition(axisR.x0 + 10, axis.y0 - qy * log10f(trh) - 15);
		treashHoldLine.text.setFillColor(sf::Color{ 150,150,150 });
		treashHoldLine.text.setFont(font);
		treashHoldLine.text.setCharacterSize(10);
		treashHoldLine.text.setString(to_string(trh));

	}


	void drawline(float ny, int lineN,sf::Color icol) {
		sf::Vertex vertexs;
		vertexs.position = sf::Vector2f(lines[lineN].x + axis.x0, -log10f(ny) * qy + axis.y0);
		//printf("\nx= %3f\ty= %3f\n", lines[lineN].x, -ny * 10 + 400.0f);
		vertexs.color = icol;

		lines[lineN].bDrawLine = true;




		if (lines[lineN].x + axis.x0 > axisR.x0) {

			for (int j = 0; j < (lines[lineN].vertex.size() - 1); j++) {
				lines[lineN].vertex[j].position.y = lines[lineN].vertex[j + 1].position.y;
			}
			lines[lineN].vertex[lines[lineN].vertex.size() - 1].position.y = vertexs.position.y;

			//cout << "\n" << lines[i].vertex.size() << endl;

		}
		else {
			lines[lineN].vertex.push_back(vertexs);
		}

		
		lines[lineN].posLine_Y[1].position.y = lines[lineN].vertex[lines[lineN].vertex.size() - 1].position.y;
		lines[lineN].posLine_Y[1].position.x = axisR.x0;
		lines[lineN].posLine_Y[1].color = sf::Color{ 50,100,50 };

		lines[lineN].posLine_Y[0].position.x = axis.x0;
		lines[lineN].posLine_Y[0].position.y = lines[lineN].posLine_Y[1].position.y;
		lines[lineN].posLine_Y[0].color = sf::Color{ 50,100,50 };



		lines[lineN].text.setPosition(axisR.x0 + 10, lines[lineN].posLine_Y[1].position.y - 15);
		lines[lineN].text.setFillColor(sf::Color{ 150,150,150 });
		lines[lineN].text.setFont(font);
		lines[lineN].text.setCharacterSize(15);
		
	}

	void initializeAxis() {
		axis.x0 = 10;
		axis.y0 = hg - 10;
		sf::Vertex dumv;

		dumv.color = sf::Color{ 100,100,100 };
		dumv.position = sf::Vector2f(0, axis.y0);
		axis.xAxis.push_back(dumv);
		dumv.position = sf::Vector2f(wg, axis.y0);
		axis.xAxis.push_back(dumv);

		dumv.color = sf::Color{ 100,100,100 };
		dumv.position = sf::Vector2f(axis.x0, 0);
		axis.yAxis.push_back(dumv);
		dumv.position = sf::Vector2f(axis.x0, hg);
		axis.yAxis.push_back(dumv);


		axisR.x0 = wg - 100;
		axisR.y0 = 0;

		dumv.color = sf::Color{ 100,100,100 };
		dumv.position = sf::Vector2f(axisR.x0, 0);
		axisR.yAxis.push_back(dumv);
		dumv.position = sf::Vector2f(axisR.x0, hg);
		axisR.yAxis.push_back(dumv);

	}

	void boundscheck() {

		bool bBoarder = false;
		int ln;
		float dx = 1;
		for (int i = 0; i < lineCount & !bBoarder; i++) {
			if (lines[i].x + axis.x0 > wg - 100) {
				bBoarder = true;
				ln = i;
				dx = (wg - axis.x0) / lines[i].x;
			}
		}

		if (bBoarder) {
			for (int i = 0; i < lineCount; i++) {
				/*
				for (int j = 0; j < lines[i].vertex.size(); j++) {
					lines[i].vertex[j].position.x = (lines[i].vertex[j].position.x-axis.x0)* dx+axis.x0;
					lines[i].x = (lines[i].vertex[j].position.x ) * dx ;
				}
				*/


				for (int j = 1; j < (lines[i].vertex.size()); j++) {
					//bufVer.position = lines[i].vertex[j].position;
					lines[i].vertex[j - 1].position = lines[i].vertex[j].position;

					//lines[i].x = (lines[i].vertex[j].position.x);
				}
				cout << "\n" << lines[i].vertex.size() << endl;
				//lines[i].vertex.pop_back();
				lines[i].x = lines[i].vertex[lines[i].vertex.size()].position.x;
			}
		}



	}


	void scale_y_max() {


		float maxY = axis.y0;

		for (int i = 0; i < lineCount; i++) {

			for (int j = 0; j < (lines[i].vertex.size()); j++) {
				if (lines[i].bdraw)
					if (lines[i].vertex[j].position.y < maxY) {
						maxY = lines[i].vertex[j].position.y;
					}
			}
		}



		if (maxY < maxH | bScale_y) {

			bScale_y = true;

			if (dqy < 1) {
				dqy += 0.01f;
			}
			else if (dqy >= 1.0f) {
				//dqy = 1;
				bScale_y = false;
				dqy = 0;
			}

			for (int i = 0; i < lineCount; i++) {

				for (int j = 0; j < lines[i].vertex.size(); j++) {
					lines[i].vertex[j].position.y = axis.y0 * (1 - (axis.y0 - maxH) / (axis.y0 - maxY))
						+ lines[i].vertex[j].position.y * (axis.y0 - maxH) / (axis.y0 - maxY);

				}


			}
			qy = qy * (axis.y0 - maxH) / (axis.y0 - maxY);

		}

	}


	void check_y() {




		float maxY = lines[0].posLine_Y[1].position.y;
		int maxI = 0;
		for (int i = 1; i < lineCount; i++) {
			if (lines[i].posLine_Y[1].position.y > maxY) {
				maxY = lines[i].posLine_Y[1].position.y;
				maxI = i;
			}

		}


		if (lines[maxI].posLine_Y[1].position.y >= axis.y0 - 10 | bCheck_y) {

			bCheck_y = true;

			if (dqy < 1) {
				dqy += 0.01f;
			}
			else if (dqy >= 1.0f) {
				//dqy = 1;
				bCheck_y = false;
				dqy = 0;
			}
			dqy = 1;
			for (int i = 0; i < lineCount; i++) {

				for (int j = 0; j < lines[i].vertex.size(); j++) {
					lines[i].vertex[j].position.y = axis.y0 * (1 - dqy * (axis.y0 - hg / 2) / (axis.y0 - maxY))
						+ lines[i].vertex[j].position.y * dqy * (axis.y0 - hg / 2) / (axis.y0 - maxY);

				}


			}
			qy = qy * dqy * (axis.y0 - hg / 2) / (axis.y0 - maxY);

		}

	}


};



