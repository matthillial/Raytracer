//UMN CSCI 5607 2D Geometry Library Homework Test Harness [HW0]
// The example file showcases some of the basic functions of the PGA library

//To compile: g++ -fsanitize=address -g -std=c++11 pga_example.cpp -o pga_example

#define _CRT_SECURE_NO_WARNINGS //Needed for MSVC to allow printf()

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>

#include "pga.h"


using namespace std;

int main(int argc, char *argv[]){
  printf("\nCSCI 5607 HW 1 - PGA Library  Examples\n");

  cout << "\nCreating Points:" << endl;
  Point2D p1 = Point2D(3,2);
  Point2D p2 = Point2D(1,1);
  Point2D p3 = Point2D(2,3);
  cout << "p1: " << string(p1) << endl;
  cout << "p2: " << string(p2) << endl;
  cout << "p3: " << string(p3) << endl;

  cout << "\nAdding and Subtraction:" << endl;
  Dir2D d1 = p2-p1;
  Point2D p4 = p1+d1;
  cout << "d1: " << string(d1) << endl;
  cout << "p4: " << string(p4) << endl;

  cout << "\nMagnitudes & Normalization:" << endl;
  Line2D l1 = Line2D(3,4,-1);
  cout << "||d1|| = " << d1.magnitude() << endl;
  cout << "d1.normalized() = " << string(d1.normalized()) << endl;
  cout << "||d1.normalized()|| = " << d1.normalized().magnitude() << endl;
  cout << "||l1|| = " << l1.magnitude() << endl;
  cout << "l1.normalized() = " << string(l1.normalized()) << endl;
  cout << "||l1.normalized()|| = " << l1.normalized().magnitude() << endl;

  cout << "\nWedge, Vee and Dot:" << endl;
  cout << "p1 v p2 : " << string(vee(p1,p2)) << endl;
  cout << "p1 v l1.normalized(): " << vee(p1,l1.normalized()) << endl;
  Line2D l2 = Line2D(5,-2,-1);
  cout << "l1 ^ l2 : " << string(wedge(l1,l2)) << endl;
  cout << "p1 . l1 : " << string(dot(p1,l1)) << endl;
  cout << "l1 . l2 : " << dot(l1,l2) << endl;
  cout << "(p1 . l1)l1 : " << string(dot(p1,l1)*l1) << endl;
  cout << "(p1 . l1)p1 : " << string(dot(p1,l1)*p1) << endl;
  cout << "p1 v p2 v p3 : " << vee(p1,p2,p3) << endl;
  cout << "sign(p1 v p2 v p3) : " << sign(vee(p1,p2,p3)) << endl;

}