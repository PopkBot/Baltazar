#include<thread>
#include<iostream>
#include <chrono>

#include <SFML/Graphics.hpp>


#include "Rocket.h"
#include "train.h"
#include "Graph.h"



using namespace std;








#define threadsCount 3 //1,2,3


//Rocket rocket;
double r = 10;

//thread tr1, tr2, tr3;

train train0;
train train1;
train train2;
train train3;

Graph graph;

bool
btrh1 = true,
btrh2 = true,
btrh3 = true;



	void thread1() {
		srand(time(0)+1);
		printf("\n=============thread1 launched==========\n");
		
		
		while (btrh1) {
			train1.control();
			train1.play_NoVisual();
			//printf("\n\ntr1 %f", train1.baltazar[0].l[4].nodes[1].bias);
		}
		printf("\n=============thread1 ceased==========\n");
	}
	
	void thread2() {
		srand(time(0)+2);
		printf("\n=============thread2 launched==========\n");
		
		//train2.initialize();
		while (btrh2) {
			train2.control();
			train2.play_NoVisual();

		}
		printf("\n=============thread2 ceased==========\n");
	}

	void thread3() {
		srand(time(0)+3);
		printf("\n=============thread3 launched==========\n");
		
		//train3.initialize();
		while (btrh3) {
			train3.control();
			train3.play_NoVisual();

		}
		printf("\n=============thread3 ceased==========\n");
	}

	void graphs() {
		graph.initialize();
	}
	







int main() {
	srand(time(0));
	int w = 1080, h = 1080;
	int mouseX = (int)(w / 2);
	int mouseY = (int)(h / 2);
	float mouseSensitivity = 3.0f;
	float speed = 0.1f;
	bool mouseHidden = false;
	bool povTower = true;
	bool keyFlag[6] = { false, false, false, false, false, false };

	float dirMoveX = 0;
	float dirMoveY = 0;
	float dirMoveZ = 0;

	sf::Vector3f pos = sf::Vector3f(-5.0f, 0.0f, 0.0f);
	sf::Clock clock;
	int framesStill = 1;


	sf::RenderWindow window(sf::VideoMode(h, w), "Tracing", sf::Style::Titlebar | sf::Style::Close); //Titlebar | sf::Style::Close
	window.setFramerateLimit(60);
	window.setPosition(sf::Vector2i(830, 0));

	sf::RenderTexture empityTexture;
	empityTexture.create(h, w);
	sf::Sprite empitySprite = sf::Sprite(empityTexture.getTexture());


	sf::Shader shader;
	shader.loadFromFile("Shader.frag", sf::Shader::Fragment);
	shader.setUniform("u_resolution", sf::Vector2f(h, w));



	train0.initialize();


	if (bdrawGraph) {
		
		graph.initializeAxis();
		graph.gridInitialize();
		graph.initializeTHLine(tresholdScore);
		train0.grLineN = 0;
		train0.grColor = sf::Color{ 220,20,60 };
		graph.lines[0].x = 0;
		graph.markerInitialize(0, sf::Color{ 220,20,60 }, "0");
		
	}







	if (Fprogramm == 1 | Fprogramm==3) {
		if (threadsCount >= 1) {
			train1.initialize();
			if (bdrawGraph) {
				train1.grLineN = 1;
				train1.grColor = sf::Color{ 255,215,0 };
				graph.lines[1].x = 0;
				graph.markerInitialize(1, sf::Color{ 255,215,0 }, "1");
			}
			thread tr1(thread1);
			tr1.detach();
		}
		if (threadsCount >= 2) {
			train2.initialize();
			if (bdrawGraph) {
				train2.grLineN = 2;
				train2.grColor = sf::Color{ 0,255,0 };
				graph.lines[2].x = 0;
				graph.markerInitialize(2, sf::Color{ 0,255,0 }, "2");
			}
			thread tr2(thread2);
			tr2.detach();
		}
		if (threadsCount >= 3) {
			train3.initialize();
			if (bdrawGraph) {
				train3.grLineN = 3;
				train3.grColor = sf::Color{ 0,0,205 };
				graph.lines[3].x = 0;
				graph.markerInitialize(3, sf::Color{ 0,0,205 }, "3");
			}
			thread tr3(thread3);
			tr3.detach();
		}

		
	}

	thread gr(graphs);
	gr.detach();



	

	

	while (window.isOpen()) {
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

		sf::Event event;

		

		while (window.pollEvent(event))
		{


			if (event.type == sf::Event::MouseWheelScrolled) {
				//std::cout « (double)event.mouseWheelScroll.delta « "\n";
				r += (double)event.mouseWheelScroll.delta;
			}

			if (event.type == sf::Event::Closed)
			{

				train0.fileOutRez(0);
				if (threadsCount >= 1) {
					train1.fileOutRez(1);
				}
				if (threadsCount >= 2) {
					train2.fileOutRez(2);
				}
				if (threadsCount >= 3) {
					train3.fileOutRez(3);
				}
				window.close();
			}
			if (event.type == sf::Event::MouseMoved)
			{

				if (mouseHidden)
				{
					int mx = event.mouseMove.x - w / 2;
					int my = event.mouseMove.y - h / 2;
					mouseX += mx;
					mouseY += my;
					sf::Mouse::setPosition(sf::Vector2i(w / 2, h / 2), window);
					if (mx != 0 || my != 0) framesStill = 1;
				}
			}
			if (event.type == sf::Event::KeyPressed) {
				if (event.key.code == sf::Keyboard::Escape) {
					window.close();
				}
				if (event.key.code == sf::Keyboard::M) {
					if (!mouseHidden) {
						framesStill = 1;
						window.setMouseCursorVisible(false);
						mouseHidden = true;
					}
					else {
						window.setMouseCursorVisible(true);
						mouseHidden = false;
					}
				}
				if (event.key.code == sf::Keyboard::N) {
					povTower = !povTower;
				}


				
				if (event.key.code == sf::Keyboard::Space) { train0.newParam(); train0.restart(); }
				if (event.key.code == sf::Keyboard::S) { train0.fileOutNN(); train0.fileOutNN_forRead(); }
				//if (event.key.code == sf::Keyboard::C) { dirMoveZ += -1; }// rocket.betta = 0;rocket.alpha = 0;rocket.gamma = 0; }
				//if (event.key.code == sf::Keyboard::W) { dirMoveX += 1; } //rocket.betta += 0.1f; }
			//	if (event.key.code == sf::Keyboard::S) { dirMoveX += -1; } //rocket.betta -= 0.1f; }
			//	if (event.key.code == sf::Keyboard::A) { dirMoveY += -1; } //rocket.gamma += 0.1f; }
			//	if (event.key.code == sf::Keyboard::D) { dirMoveY += 1; } //rocket.gamma -= 0.1f;}
				//if (event.key.code == sf::Keyboard::Q) { rocket.bTurnX = 1; } //rocket.alpha += 0.1f; }
				//if (event.key.code == sf::Keyboard::E) { rocket.bTurnX = -1; } //rocket.alpha -= 0.1f; }
				//printf("bthrust =%d\n", rocket.bThrust);
				
			}
			else if (event.type == sf::Event::KeyReleased)
			{
				
				//if (event.key.code == sf::Keyboard::Space)dirMoveZ = 0;
				//if (event.key.code == sf::Keyboard::C)dirMoveZ = 0; 
				//if (event.key.code == sf::Keyboard::W)dirMoveX = 0;
				//if (event.key.code == sf::Keyboard::S)dirMoveX = 0;
				//if (event.key.code == sf::Keyboard::A)dirMoveY = 0;
				//if (event.key.code == sf::Keyboard::D)dirMoveY = 0;
				//if (event.key.code == sf::Keyboard::Q)rocket.bTurnX = 0;
				//if (event.key.code == sf::Keyboard::E)rocket.bTurnX = 0;

				if (event.key.code == sf::Keyboard::Num1) { train0.fileOutNN(); train0.fileOutNN_forRead(); printf("\n train 0 SAVED\n"); }
				if (event.key.code == sf::Keyboard::Num2) { train1.fileOutNN(); train1.fileOutNN_forRead(); printf("\n train 1 SAVED\n"); }
				if (event.key.code == sf::Keyboard::Num3) { train2.fileOutNN(); train2.fileOutNN_forRead(); printf("\n train 2 SAVED\n"); }
				if (event.key.code == sf::Keyboard::Num4) { train3.fileOutNN(); train3.fileOutNN_forRead(); printf("\n train 3 SAVED\n"); }

			}
			
		}
		float mx = ((float)mouseX / w - 0.5f) * mouseSensitivity;
		float my = ((float)mouseY / h - 0.5f) * mouseSensitivity;
		if (mouseHidden) {


			sf::Vector3f dir = sf::Vector3f(0.0f, 0.0f, 0.0f);
			sf::Vector3f dirTemp;
			

			dirTemp.z += -(r+dirMoveZ) * sin(my);
			dirTemp.x += -(r+dirMoveX) * cos(my) * cos(mx);
			dirTemp.y += -(r+dirMoveY) * cos(my) * sin(mx);
			if (!keyFlag[0]) {
				dir.x = dirTemp.x;
				dir.y = dirTemp.y;
				dir.z = dirTemp.z;
			}
			else {
				dir.x = dirTemp.x * cos(mx) - dirTemp.y * sin(mx);
				dir.y = dirTemp.x * sin(mx) + dirTemp.y * cos(mx);
				dir.z = dirTemp.z;

			}
			pos = dir;
		}
		train0.control();
		train0.play(shader);
	

		if (povTower) shader.setUniform("u_pos", pos + sf::Vector3f{0.0f,0.0f,-120.0f});
		else shader.setUniform("u_pos", pos + train0.rocket[0].rOA);
		


		shader.setUniform("u_mouse", sf::Vector2f(mx, my));
		shader.setUniform("u_time", clock.getElapsedTime().asSeconds());
		
		window.draw(empitySprite, &shader);

		window.display();
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		//std::cout « (float)round(round((end - begin).count() /1000)/1000)/1000 « "\n";
		//rocket.dt = (float)round(round((end - begin).count() / 1000) / 1000) / 1000;
		//std::cout « "Time difference = " « std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() « "[µs]" « std::endl;

	}

	btrh1 = false;
	btrh2 = false;
	btrh3 = false;

	//train0.fileOutNN();
	//train0.fileOutNN_forRead();
	


	return 0;
}