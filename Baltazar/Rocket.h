#pragma once

#include <SFML/Graphics.hpp>
#include<iostream>


//#include"Vec3.h"



#define kRock 0.2f
#define wind10 0.5f
#define windfactor 0.25f
#define hr 2


class Rocket
{
public:

	const double kalpha = round(rand() % 360) / 3.1415;

	sf::Vector3f rOC, rAB = { hr,0,0 }, rOM, rOA = { 1,1,-100 }, rOB, rAC, rVecCtrl1 = { 0,1,0 }, rVecCtrl2 = { 0,0,1 }, vecThrust, track[100], vecX = { 1,0,0 }, vecY = { 0,1,0 }, vecZ = { 0,0,1 };
	float
		ra = 0.5f,
		rb = 0.2f,
		h = hr,
		h2 = 0.5,
		lctrl = 0.2f,
		lvecCtlr = 0.5f;


	int
		trackCount,
		engIgn=1;

	float fireDraw[7]= { 0,0,0,0 ,0,0,0 };
	float fire[7] = { 0,0,0,0 ,0,0,0 };
	sf::Vector3f vecCtrl[7];



	float
		dt = 0.016f,
		t = 0,
		m0,
		mf,
		dm,
		m = 1,
		jox = 1,
		joy = 1,
		joz = 1,
		dax,
		day,
		daz,
		dEox,
		dEoy,
		dEoz,
		alpha ,
		betta ,
		gamma ,
		bThrust = 0,
		
		phi = 0,
		tetta = 0,
		fctrl[7],
		//vx = 0,
		//vy,
		//vz,
		//wx = 0,
		//wy = 0,
		//wz = 0,
		//tralpha,
		//trbetta,
		//xm = 0,
		//ym = 0,
		//zm = 0,
		thrust = 20,
		ct = 1.0f,
		cx = 0.2f,
		cy = 0.5f,
		cz = 0.5f,
		sX = 0.01f,
		sY = 0.01f,
		sZ = 0.01f;


	sf::Vector3f

		g = { -9.81f,0,0 },
		fvec[7],
		fA,
		vm ,
		rm,
		wm ,
		vA ,
		vCD = { h * 0.25f,0,0 }, //vектор от ЦМ к ЦД
		da ;

	int
		
		bStraifY = 0,
		bStraifZ = 0,
		bTurnX = 0,
		bTurnY = 0,
		bTurnZ = 0;


	bool dead,
		engRun=false
		;

	double
		cfuel = 0,
		tfuel=0;




	void draw(sf::Shader& shader) { //sf::Shader& shader,sf::Texture textureGrass


	//rAB = coordTransorm(vecFullUnTurn({h,0,0},alpha,betta,gamma));

		rOB = rOA + coordTransorm(rAB);




		vecCtrl[1] = coordTransorm(normalize(-cross(vecZ, rAB)) * lctrl + rAB) + rOA;
		vecCtrl[4] = coordTransorm(normalize(-cross(rAB, vecZ)) * lctrl + rAB) + rOA;
		vecCtrl[2] = coordTransorm(normalize(cross(vecY, rAB)) * lvecCtlr) + vecCtrl[1];
		vecCtrl[3] = coordTransorm(normalize(cross(rAB, vecY)) * lvecCtlr) + vecCtrl[1];
		vecCtrl[5] = coordTransorm(normalize(cross(vecY, rAB)) * lvecCtlr) + vecCtrl[4];
		vecCtrl[6] = coordTransorm(normalize(cross(rAB, vecY)) * lvecCtlr) + vecCtrl[4];
		vecThrust = coordTransorm(normalize(-vecX * cos(tetta) * cos(phi) + vecY * cos(tetta) * sin(phi) + vecZ * sin(tetta)) * thrust / 20.0f) + rOA;
		//vecThrust = coordTransorm(normalize({ - cos(tetta) * cos(phi), cos(tetta) * sin(phi), sin(tetta) })*thrust/20.0f) + rOA;
		//{ -thrust / 20 * cos(tetta) * cos(phi), thrust / 20 * cos(tetta) * sin(phi), thrust / 20 * sin(tetta) }
		/*
		vecCtrl[1] = coordTransorm(vecFullUnTurn((rVecCtrl1 * lctrl), alpha, betta, gamma)) + rAB + rOA;
		vecCtrl[2] = coordTransorm(vecFullUnTurn((rVecCtrl2 * lvecCtlr + rVecCtrl1 * lctrl), alpha, betta, gamma)) + rAB + rOA;
		vecCtrl[3] = coordTransorm(vecFullUnTurn((-rVecCtrl2 * lvecCtlr + rVecCtrl1 * lctrl), alpha, betta, gamma)) + rAB + rOA;
		vecCtrl[4] = coordTransorm(vecFullUnTurn((-rVecCtrl1 * lctrl), alpha, betta, gamma)) + rAB + rOA;
		vecCtrl[5] = coordTransorm(vecFullUnTurn((rVecCtrl2 * lvecCtlr - rVecCtrl1 * lctrl), alpha, betta, gamma)) + rAB + rOA;
		vecCtrl[6] = coordTransorm(vecFullUnTurn((-rVecCtrl2 * lvecCtlr - rVecCtrl1 * lctrl), alpha, betta, gamma)) + rAB + rOA;

		vecThrust = coordTransorm(vecFullUnTurn({ -thrust/20 * cos(tetta) * cos(phi),thrust/20 * cos(tetta) * sin(phi),thrust/20 * sin(tetta) }, alpha, betta, gamma)) + rOA;;
		*/


		//shader.setUniform("u_textureGrass", textureGrass);
		shader.setUniform("u_thrust", vecThrust);
		shader.setUniform("u_coneROA", rOA);
		shader.setUniform("u_coneROB", rOB);
		shader.setUniform("u_lvecCtlr", lvecCtlr);
		shader.setUniformArray("u_vecCtrl", vecCtrl, 7);
		shader.setUniformArray("u_fire", fireDraw, 7);
		shader.setUniform("u_ra", ra);
		shader.setUniform("u_rb", rb);


	}

	void preDraw() {
		rOB = rOA + coordTransorm(rAB);
		vecCtrl[1] = coordTransorm(normalize(-cross(vecZ, rAB)) * lctrl + rAB) + rOA;
		vecCtrl[4] = coordTransorm(normalize(-cross(rAB, vecZ)) * lctrl + rAB) + rOA;
		vecCtrl[2] = coordTransorm(normalize(cross(vecY, rAB)) * lvecCtlr) + vecCtrl[1];
		vecCtrl[3] = coordTransorm(normalize(cross(rAB, vecY)) * lvecCtlr) + vecCtrl[1];
		vecCtrl[5] = coordTransorm(normalize(cross(vecY, rAB)) * lvecCtlr) + vecCtrl[4];
		vecCtrl[6] = coordTransorm(normalize(cross(rAB, vecY)) * lvecCtlr) + vecCtrl[4];
		vecThrust = coordTransorm(normalize(-vecX * cos(tetta) * cos(phi) + vecY * cos(tetta) * sin(phi) + vecZ * sin(tetta)) * thrust / 20.0f) + rOA;
	}

	/*
	void spinning() {
		alpha = 0.1f;
		//betta =0.79f;
		//gamma +=0.001f;

		/*if (rOA.x > 20) rOA.x = 0;
		rOA.x += 0.05;
		rOA.y = rOA.x * rOA.x / 20;
		rOA.z = (-rOA.x * 2);
		*/
		/*if (rOA.x * rOA.x + rOA.y * rOA.y + rOA.z * rOA.z < 1000000)
		rOA += coordTransorm(vecFullUnTurn(sf::Vector3f(1, 0,0), alpha,

		betta, gamma));
		else rOA = sf::Vector3f(0, 0, 0);
		
		printf("%5f\t%5f\t%5f\t%5f\n", rAB.x, rAB.y, rAB.z, sqrtf(pow(rAB.x, 2) + pow(rAB.y, 2) + pow(rAB.z, 2)));
		rAB = qatTurn(rAB, normalize({ 1,1,1 }), alpha);
		vecZ = qatTurn(vecZ, normalize({ 1,1,1 }), alpha);
		vecY = qatTurn(vecY, normalize({ 1,1,1 }), alpha);
		vecX = qatTurn(vecX, normalize({ 1,1,1 }), alpha);
	}
	*/

	void initializeTurn() {
		vecX = { 1,0,0 };
		vecY = { 0,1,0 };
		vecZ = { 0,0,1 };
		rAB = { hr,0,0 };
		vecZ = qatTurn(vecZ, normalize(vecX), alpha);
		vecY = qatTurn(vecY, normalize(vecX), alpha);
		rAB = qatTurn(rAB, normalize(vecX), alpha);

		vecZ = qatTurn(vecZ, normalize(vecY), betta);
		vecX = qatTurn(vecX, normalize(vecY), betta);
		rAB = qatTurn(rAB, normalize(vecY), betta);

		vecX = qatTurn(vecX, normalize(vecZ), gamma);
		vecY = qatTurn(vecY, normalize(vecZ),gamma);
		rAB = qatTurn(rAB, normalize(vecZ),gamma);

	}

	void dynamic() {

		vecX = normalize(vecX);
		vecY = normalize(vecY);
		vecZ = normalize(vecZ);
		//vm = vecFullTurn(vm, alpha, betta, gamma);

		vA.y = wind10 * sin(kalpha) * pow(abs(rOA.z) / 10, kRock) + windfactor * round(rand() % 20 - 10) / 10;
		vA.z = wind10 * cos(kalpha) * pow(abs(rOA.z) / 10, kRock) + windfactor * round(rand() % 20 - 10) / 10;
		vA.x = windfactor * round(rand() % 4 - 2) / 10.0f;



		fA.x = cx * 0.635f * abs(dot(vA - vm, vecX)) * dot(vA - vm, vecX) * sX;
		fA.y = cy * 0.635f * abs(dot(vA - vm, vecY)) * dot(vA - vm, vecY) * sY;
		fA.z = cz * 0.635f * abs(dot(vA - vm, vecZ)) * dot(vA - vm, vecZ) * sZ;
		fA = fA.x * vecX + fA.y * vecY + fA.z * vecZ;

		fA = { 0,0,0 };

		/*
		dax = vecFullTurn(g, alpha, betta, gamma).x+(thrust*fire[0]*cos(-tetta)*cos(phi) +fA.x)/m;
		day = vecFullTurn(g, alpha, betta, gamma).y+(thrust*fire[0]*cos(-tetta)*sin(-phi) +fA.y)/m; //+ct*0*(-fire[1]+fire[4])
		daz = vecFullTurn(g, alpha, betta, gamma).z+(thrust*fire[0]*sin(-tetta) +fA.z)/m; //+ct*0*(-fire[2]-fire[5]+fire[3]+fire[6])
		vm.x += dax * dt;
		vm.y += day * dt;
		vm.z += daz * dt;
		*/
		/*
		dax =(thrust * fire[0] * cos(-tetta) * cos(phi)) / m;
		day =(thrust * fire[0] * cos(-tetta) * sin(-phi)) / m; //+ct*0*(-fire[1]+fire[4])
		daz = (thrust * fire[0] * sin(-tetta)) / m; //+ct*0*(-fire[2]-fire[5]+fire[3]+fire[6])
		*/

		sf::Vector3f da1;
		da1 = normalize(vecX * cos(tetta) * cos(phi) - vecY * cos(tetta) * sin(phi) - vecZ * sin(tetta)) * thrust * fire[0] + fA;
		da = sf::Vector3f{ dot(da1,{1,0,0}), dot(da1,{0,1,0}), dot(da1,{0,0,1}) } + g;
		vm += da * dt;

		rOA += coordTransorm(vm * dt);


		dEox = (ct * (-fire[2] - fire[6] + fire[3] + fire[5]) * lctrl + dot(fA, vecY) * vCD.z - dot(fA, vecZ) * vCD.y) / jox;//+ fA.y*vCD.z-fA.z*vCD.y
		wm.x += dEox * dt ;
		vecZ = qatTurn(vecZ, normalize(vecX), wm.x * dt);
		vecY = qatTurn(vecY, normalize(vecX), wm.x * dt);
		rAB = qatTurn(rAB, normalize(vecX), wm.x * dt);

		dEoy = (ct * (fire[3] + fire[6] - fire[2] - fire[5] - dot(fA, vecZ) * vCD.x + dot(fA, vecX) * vCD.z) * h) / joy;//+ fA.z*vCD.x-fA.x*vCD.z
		wm.y += dEoy * dt ;
		vecZ = qatTurn(vecZ, normalize(vecY), wm.y * dt);
		vecX = qatTurn(vecX, normalize(vecY), wm.y * dt);
		rAB = qatTurn(rAB, normalize(vecY), wm.y * dt);

		dEoz = (ct * (-fire[4] + fire[1]) * h - dot(fA, vecX) * vCD.y + dot(fA, vecY) * vCD.x) / joz;//+ fA.x*vCD.y-fA.y*vCD.x
		wm.z += dEoz * dt ;
		vecX = qatTurn(vecX, normalize(vecZ), wm.z * dt);
		vecY = qatTurn(vecY, normalize(vecZ), wm.z * dt);
		rAB = qatTurn(rAB, normalize(vecZ), wm.z * dt);

		alpha += wm.x * dt;
		betta += wm.y * dt;
		gamma += wm.z * dt;


		/*
		sf::Vector3f lvec = vecFullUnTurn({ 0,dEoy,dEoz },alpha,betta,gamma) ; // по часовой + против часовой -

		wm.x += dEox*dt;
		wm.y += lvec.y*dt ;
		wm.z += lvec.z*dt ;

		alpha += wm.x * dt;
		betta += wm.y * dt;
		gamma += wm.z * dt;
		*/

		t += dt;

		//printf("%10f\t%1f\t%1f\t%3f\n", dEox, dEoy, dEoz, -rOA.z);
		//if(rOA.z<0 ) printf("%5f %5f \n",t,-rOA.z);
	}




	void control() {

		if (!engRun & bThrust > 0 & engIgn > 0) {
			engIgn -= 1;
			engRun = true;
		}
		if (engRun & bThrust <= 0) {
			engRun = false;
		}
		if (engIgn < 0) engIgn = 0;



		if (engRun) {
			fire[0] = bThrust;
			fireDraw[0] = 1;
		}
		else {
			fire[0] = 0;
			fireDraw[0] = 0;
		}


		fire[1] = -bTurnZ;
		if (fire[1] < 0) fire[1] = 0; 


		fire[2] = +bTurnY + bTurnX;
		if (fire[2] < 0)fire[2] = 0;

		fire[3] = -bTurnY - bTurnX;
		if (fire[3] < 0)fire[3] = 0;


		fire[4] = +bTurnZ;
		if (fire[4] < 0) fire[4] = 0;

		fire[5] = +bTurnY - bTurnX;
		if (fire[5] < 0)fire[5] = 0;

		fire[6] = -bTurnY + bTurnX;
		if (fire[6] < 0)fire[6] = 0;

		for (int i = 1; i < 7; i++) {
			fireDraw[i] = fire[i];
		}

	}



	void traektory(sf::Shader& shader) {
		float ft = falltime(rOA.z);


		if (ft > 0) {
			trackCount = abs((int)round(ft) / 2) + 20;

			track[0].x = rOA.x - vm.y * ft;
			track[0].y = rOA.y - vm.z * ft;
			track[0].z = 0;
			if (trackCount > 98)trackCount = 98;
			for (int i = 1; i < trackCount; i++) {


				track[i].x = rOA.x - vm.y * ft / (trackCount)*i;
				track[i].y = rOA.y - vm.z * ft / (trackCount)*i;
				track[i].z = rOA.z - vm.x * ft / (trackCount)*i - g.x * powf(ft / (trackCount)*i, 2) / 2;

			}
			shader.setUniformArray("u_track", track, trackCount);
			shader.setUniform("u_trackCount", (float)trackCount);
			shader.setUniform("u_Avec", coordTransorm(sf::Vector3f{ dot(fA,{1,0,0}), dot(fA,{0,1,0}), dot(fA,{0,0,1}) })); //vecFullUnTurn(fA, alpha, betta, gamma)	sf::Vector3f{ dot(fA,{1,0,0}), dot(fA,{0,1,0}), dot(fA,{0,0,1}) })
			shader.setUniform("u_CDvec", coordTransorm(normalize(rAB) * vecLength(vCD)));//vecFullUnTurn(vCD,alpha,betta,gamma)
		}
	}


	float falltime(float z0) {
		float ft;
		if (z0 * 0 - rOA.z > 0)
			ft = (vm.x + sqrtf(vm.x * vm.x + 2 * g.x * (-z0 * 0 + rOA.z))) / (-g.x);
		else
			ft = 0;

		return ft;
	}
	sf::Vector3f xUnTurn(sf::Vector3f vec, double alpha0) {
		sf::Vector3f newVec;
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = x0;
		newVec.y = y0 * cos(alpha0) + z0 * sin(alpha0);
		newVec.z = z0 * cos(alpha0) - y0 * sin(alpha0);
		return newVec;
	}
	sf::Vector3f yUnTurn(sf::Vector3f vec, double alpha0) {
		sf::Vector3f newVec;
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = x0 * cos(alpha0) - z0 * sin(alpha0);
		newVec.y = y0;
		newVec.z = z0 * cos(alpha0) + x0 * sin(alpha0);
		return newVec;
	}
	sf::Vector3f zUnTurn(sf::Vector3f vec, double alpha0) {
		sf::Vector3f newVec;
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = x0 * cos(alpha0) + y0 * sin(alpha0);
		newVec.y = y0 * cos(alpha0) - x0 * sin(alpha0);
		newVec.z = z0;
		return newVec;
	}
	sf::Vector3f xTurn(sf::Vector3f vec, double alpha0) {
		sf::Vector3f newVec;
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = x0;
		newVec.y = y0 * cos(alpha0) - z0 * sin(alpha0);
		newVec.z = z0 * cos(alpha0) + y0 * sin(alpha0);
		return newVec;
	}
	sf::Vector3f yTurn(sf::Vector3f vec, double alpha0) {
		sf::Vector3f newVec;
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = x0 * cos(alpha0) + z0 * sin(alpha0);
		newVec.y = y0;
		newVec.z = z0 * cos(alpha0) - x0 * sin(alpha0);
		return newVec;
	}
	sf::Vector3f zTurn(sf::Vector3f vec, double alpha0) {
		sf::Vector3f newVec;
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = x0 * cos(alpha0) - y0 * sin(alpha0);
		newVec.y = y0 * cos(alpha0) + x0 * sin(alpha0);
		newVec.z = z0;
		return newVec;
	}
	sf::Vector3f fullUnTurn(double x0, double y0, double z0, double a0, double b0, double g0) {
		sf::Vector3f newVec;
		newVec.x = x0 * cos(b0) * cos(g0) - z0 * sin(b0) + y0 * cos(b0) * sin(g0);
		newVec.y = y0 * (cos(a0) * cos(g0) + sin(a0) * sin(b0) * sin(g0)) - x0 * (cos(a0) * sin(g0) - cos(g0) * sin(a0) * sin(b0)) + z0 * cos(b0) * sin(a0);
		newVec.z = x0 * (sin(a0) * sin(g0) + cos(a0) * cos(g0) * sin(b0)) - y0 * (cos(g0) * sin(a0) - cos(a0) * sin(b0) * sin(g0)) + z0 * cos(a0) * cos(b0);
		return newVec;
	}
	sf::Vector3f vecFullUnTurn(sf::Vector3f vec, double a0, double b0, double g0) {
		sf::Vector3f newVec = vec;

		double x0 = vec.x, y0 = vec.y, z0 = vec.z;

		//newVec.x = z0 * (sin(a0) * sin(g0) + cos(a0) * sin(b0) * cos(g0)) - y0 * (sin(a0) * cos(g0) - cos(a0) * sin(b0) * sin(g0)) + x0 * cos(a0) * cos(b0);
		//newVec.y = y0 * (cos(a0) * cos(g0) + sin(a0) * sin(b0) * sin(g0)) - z0 * (cos(a0) * sin(g0) - cos(g0) * sin(a0) * sin(b0)) + x0 * cos(b0) * sin(a0);
		//newVec.z = y0 * cos(b0) * sin(g0) - x0 * sin(b0) + z0 * cos(b0) * cos(g0);

		newVec = zUnTurn(yUnTurn(xUnTurn(newVec, a0), b0), g0);

		return newVec;
	}
	sf::Vector3f vecFullTurn(sf::Vector3f vec, double a0, double b0, double g0) {
		sf::Vector3f newVec = vec;
		/*
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = z0

		* sin(b0) + x0 * cos(b0) * cos(g0) - y0 * cos(b0) * sin(g0);
		newVec.y = x0 * (cos(a0) * sin(g0) + sin(a0) * sin(b0) * cos(g0)) + y0 * (cos(a0) * cos(g0) - sin(a0) * sin(b0) * sin(g0)) - z0 * sin(a0) * cos(b0);
		newVec.z = x0*(sin(a0)*sin(g0)-cos(a0)*sin(b0)*cos(g0))+y0*(sin(a0)*cos(g0)+cos(a0)*sin(b0)*sin(g0))+z0*cos(a0)*cos(b0);
		*/
		//newVec = xTurn(yTurn(zTurn(newVec, g0),b0), a0);
		newVec = xTurn(yTurn(zTurn(newVec, g0), b0), a0);
		return newVec;
	}
	sf::Vector3f coordTransorm(sf::Vector3f vec) {
		sf::Vector3f newVec;
		newVec.x = -vec.y;
		newVec.y = -vec.z;
		newVec.z = -vec.x;
		return newVec;
	}



	sf::Vector3f qatTurn(sf::Vector3f p1, sf::Vector3f u, double a) {


		sf::Vector3f p = p1;

		double t = p.x * u.x + p.y * u.y + p.z * u.z;

		p.x = p.x * cos(a) + (p.z * u.y - u.z * p.y) * sin(a / 2) * cos(a / 2) + p.x * t * pow(sin(a / 2), 2);
		p.y = p.y * cos(a) + (p.x * u.z - u.x * p.z) * sin(a / 2) * cos(a / 2) + p.y * t * pow(sin(a / 2), 2);
		p.z = p.z * cos(a) + (p.y * u.x - u.y * p.x) * sin(a / 2) * cos(a / 2) + p.z * t * pow(sin(a / 2), 2);
		p = normalize(p) * vecLength(p1);
		return p;
	}

	sf::Vector3f normalize(sf::Vector3f vec) {

		sf::Vector3f newVec = vec;

		newVec = newVec / vecLength(newVec);
		return newVec;
	}

	sf::Vector3f vecRecoverLength(sf::Vector3f vec, float length) {
		sf::Vector3f newVec = vec;
		newVec = normalize(newVec);
		newVec = newVec * length;
		return newVec;

	}
	float vecLength(sf::Vector3f vec) {
		float length = sqrtf(pow(vec.x, 2) + pow(vec.y, 2) + pow(vec.z, 2));
		return length;
	}

	sf::Vector3f cross(sf::Vector3f vec1, sf::Vector3f vec2) {
		sf::Vector3f nv1 = vec1, nv2 = vec2, newVec;

		newVec.x = nv1.y * nv2.z - nv1.z * nv2.y;
		newVec.y = nv1.z * nv2.x - nv1.x * nv2.z;
		newVec.z = nv1.x * nv2.y - nv1.y * nv2.x;
		return newVec;
	}

	float dot(sf::Vector3f vec1, sf::Vector3f vec2) {
		float d = vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
		return d;
	}


};