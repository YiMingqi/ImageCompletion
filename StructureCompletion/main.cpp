// StructurePropagation.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include "StructurePropagation.h"
#include "OpenCvUtility.h"
#include <fstream>
Point points[2]={Point(-1,-1),Point(-1,-1)};
int points_i=0;
Mat3b img;
Mat1b mask;
Mat3b result;
Mat3b result_copy;
vector<vector<Point>> PointsList;
Point prev_pt(-1,-1);
vector<vector<Point>> mousepoints;
vector<Point> curvePoints;
int blocksize=20;
int samplestep=2;
bool iscurve=true;
ofstream file;
void onmouse(int event,int x,int y,int flags,void* parm)
{
	if(!iscurve)
	{
		if (event!=CV_EVENT_LBUTTONDOWN)
			return;
		points[points_i].x=x;
		points[points_i].y=y;
		//cout<<x<<"  "<<y<<"  "<<points[0]<<" "<<points[1]<<endl;
		points_i=(points_i+1)%2;
		if (points[0].x!=-1&&points[1].x!=-1&&points_i==0)
		{
			// result.copyTo(result_copy);
			file<<points[0]<<"  "<<points[1]<<endl;
			vector<Point> line;
			LineInterpolation(points, line);
			PointsList.push_back(line);
			DrawPoints(line, result_copy, Scalar(255, 0, 255), 1);//×ÏÉ«
			circle(result_copy,points[0],3,Scalar(255,0,0),CV_FILLED);//À¶É«
			circle(result_copy,points[1],3,Scalar(255,0,0),CV_FILLED);
			
			// rectangle(result_copy,RectByCenter(PointsList[0],blocksize),CV_RGB(255,0,0),2);
			imshow("img",result_copy);
		}
	}
	else
	{
		if (event == CV_EVENT_LBUTTONUP) {
			prev_pt = cvPoint(-1, -1);
			mousepoints.push_back(curvePoints);
			curvePoints = vector<Point>();
		}
		if( event == CV_EVENT_LBUTTONDOWN )
		{
			prev_pt = cvPoint(x,y);
			// rectangle(result_copy,RectByCenter(prev_pt,blocksize),CV_RGB(255,0,0),2);
			curvePoints.push_back(prev_pt);
		}
		else if( event == CV_EVENT_MOUSEMOVE && (flags & CV_EVENT_FLAG_LBUTTON) )
		{
			CvPoint pt = cvPoint(x,y);
			curvePoints.push_back(pt);
			if( prev_pt.x < 0 )
				prev_pt = pt;
			//cvLine( inpaint_mask, prev_pt, pt, cvScalarAll(255), 5, 8, 0 );
			line( result_copy, prev_pt, pt, cvScalarAll(255), 1, 8, 0 );
			prev_pt = pt;
			imshow( "img", result_copy );
		}
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	// img=imread("img.jpg",1);
	// mask=imread("mask.bmp",0);

	img=imread("K.png",1);
	mask=imread("K.bmp",0);

	threshold(mask,mask,125,255,CV_THRESH_BINARY_INV);
	result.zeros(img.size());
	img.copyTo(result,mask);
	namedWindow("img");
	//namedWindow("mask");
	createTrackbar("BlockSize","img",&blocksize,50);
	createTrackbar("SampleStep","img",&samplestep,20);
	int iscurve_temp=true;
	createTrackbar("iscurve","img",&iscurve_temp, 1);
	setMouseCallback("img",onmouse);
	imshow("img",result);
	//imshow("mask",mask);
	StructurePropagation SP;
	SP.SetParm(blocksize,samplestep,iscurve);
	Mat3b Local_Result_Copy(result.size());
	result.copyTo(Local_Result_Copy);
	result.copyTo(result_copy);
	//file.open("test.txt");

	for (;;)
	{
		iscurve=iscurve_temp;
		char c=cvWaitKey(10);
		if (c==27)
			break;
		else if (c=='s')
		{
			file<<blocksize<<"   "<<samplestep<<endl;
			if (iscurve) {
				PointsList.resize(mousepoints.size());
				for (int i = 0; i < mousepoints.size(); i++) {
					Wang_GetCurve(mousepoints[i], PointsList[i]);
				}
			}

			for (int i = 0; i < PointsList.size(); i++) {
				DrawPoints(PointsList[i], img, CV_RGB(255, 0, 0), 1);
			}
			SP.SetParm(blocksize,samplestep,iscurve);
			SP.Run(mask,result,PointsList,Local_Result_Copy);
			imshow("img", Local_Result_Copy);
			Local_Result_Copy.copyTo(result_copy);
			PointsList.clear();
			mousepoints.clear();
		}
		else if (c=='r')
		{
			// clear lines
			result.copyTo(result_copy);
			result.copyTo(Local_Result_Copy);
			PointsList.clear();
			mousepoints.clear();
			imshow("img",result_copy);
		}
		else if (c=='a')
		{
			// save results
			imwrite("result.jpg",Local_Result_Copy);
		}
		else if(c=='e')
		{
			PointsList.resize(mousepoints.size());
			for (int i = 0; i < mousepoints.size(); i++) {
				Wang_GetCurve(mousepoints[i], PointsList[i]);
				DrawPoints(PointsList[i], result_copy, CV_RGB(255, 0, 0), 1);
			}
			imshow("img",result_copy);
		}
	}
	file.close();
	return 0;
}

