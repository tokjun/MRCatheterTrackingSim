/*=========================================================================

  Program:   OpenIGTLink -- Example for Tracker Client Program
  Language:  C++

  Copyright (c) Insight Software Consortium. All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <iomanip>

#include "igtlOSUtil.h"
#include "igtlTrackingDataMessage.h"
#include "igtlClientSocket.h"
#include "igtlMath.h"


typedef struct {
  double ts;
  double tip1X;
  double tip1Y;
  double tip1Z;
  double tip1NX;
  double tip1NY;
  double tip1NZ;
  double coil1X;
  double coil1Y;
  double coil1Z;
  double coil2X;
  double coil2Y;
  double coil2Z;
  double tip2X;
  double tip2Y;
  double tip2Z;
  double tip2NX;
  double tip2NY;
  double tip2NZ;
  double coil3X;
  double coil3Y;
  double coil3Z;
  double coil4X;
  double coil4Y;
  double coil4Z;
} TrackingData;

typedef struct {
  double ts;
  igtl::Matrix4x4 matrixCoil1;
  igtl::Matrix4x4 matrixCoil2;
  igtl::Matrix4x4 matrixCoil3;
  igtl::Matrix4x4 matrixCoil4;
} MatrixData;

typedef std::vector<TrackingData> TrackingDataList;
typedef std::vector<MatrixData> MatrixList;

void  ReadFile(std::string filename, TrackingDataList& coordinates);
void  ConvertTrackingData(TrackingDataList& coordinates, MatrixList& matrices);
int   SendTrackingData(igtl::ClientSocket::Pointer& socket, igtl::TrackingDataMessage::Pointer& trackingMsg, MatrixData& mat);
int   SendTrackingData(igtl::ClientSocket::Pointer& socket, igtl::TrackingDataMessage::Pointer& trackingMsg);
void  GetRandomTestMatrix(igtl::Matrix4x4& matrix, float phi, float theta);

int main(int argc, char* argv[])
{
  //------------------------------------------------------------
  // Parse Arguments

  if (argc != 5) // check number of arguments
    {
    // If not correct, print usage
    std::cerr << "Usage: " << argv[0] << " <hostname> <port> <fps>"    << std::endl;
    std::cerr << "    <hostname> : IP or host name"                    << std::endl;
    std::cerr << "    <port>     : Port # (18944 in Slicer default)"   << std::endl;
    std::cerr << "    <fps>      : Frequency (fps) to send coordinate" << std::endl;
    std::cerr << "    <file>     : Tracking file" << std::endl;
    exit(0);
    }

  char*  hostname = argv[1];
  int    port     = atoi(argv[2]);
  double fps      = atof(argv[3]);
  int    interval = (int) (1000.0 / fps);
  std::string filename = argv[4];


  //------------------------------------------------------------
  // Load Tracking Data
  
  TrackingDataList coordinates;
  MatrixList matrices;

  ReadFile(filename, coordinates);
  ConvertTrackingData(coordinates, matrices);
  
  //------------------------------------------------------------
  // Establish Connection

  igtl::ClientSocket::Pointer socket;
  socket = igtl::ClientSocket::New();
  int r = socket->ConnectToServer(hostname, port);

  if (r != 0)
    {
    std::cerr << "Cannot connect to the server." << std::endl;
    exit(0);
    }

  //------------------------------------------------------------
  // Allocate TrackingData Message Class
  //
  // NOTE: TrackingDataElement class instances are allocated
  //       before the loop starts to avoid reallocation
  //       in each image transfer.
  
  igtl::TrackingDataMessage::Pointer trackingMsg;
  trackingMsg = igtl::TrackingDataMessage::New();

  igtl::TrackingDataElement::Pointer trackElement[4];
  
  for (int coil = 0; coil < 4; coil ++)
    {
    std::stringstream ss;
    ss << "Tracker" << std::setfill('0') << std::setw(2) << coil;
    trackElement[coil] = igtl::TrackingDataElement::New();
    trackElement[coil]->SetName(ss.str().c_str());
    trackElement[coil]->SetType(igtl::TrackingDataElement::TYPE_3D);
    trackingMsg->AddTrackingDataElement(trackElement[coil]);
    }
  
  //------------------------------------------------------------
  // Loop
  MatrixList::iterator iter;
  for (iter = matrices.begin(); iter != matrices.end(); iter ++)
    {
    trackingMsg->SetDeviceName("Tracker");
    //SendTrackingData(socket, trackingMsg);
    SendTrackingData(socket, trackingMsg, *iter);
    igtl::Sleep(interval);
    }
      
}


void ReadFile(std::string filename, TrackingDataList& coordinates)
{
  std::ifstream ifs;

  ifs.open(filename);
  if (!ifs.is_open())
    {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(0);
    }
  
  while (!ifs.eof())
    {
    TrackingData pt;

    if (ifs
        >> pt.ts
        >> pt.tip1X  >> pt.tip1Y  >> pt.tip1Z
        >> pt.tip1NX >> pt.tip1NY >> pt.tip1NZ
        >> pt.coil1X >> pt.coil1Y >> pt.coil1Z
        >> pt.coil2X >> pt.coil2Y >> pt.coil2Z
        >> pt.tip2X  >> pt.tip2Y  >> pt.tip2Z
        >> pt.tip2NX >> pt.tip2NY >> pt.tip2NZ
        >> pt.coil3X >> pt.coil3Y >> pt.coil3Z
        >> pt.coil4X >> pt.coil4Y >> pt.coil4Z)
      {
      // TODO: Need LPS to RAS conversion?
      
      coordinates.push_back(pt);
      }
    else
      {
      std::cerr << "Read: " << coordinates.size() << " coordinates from file." << std::endl;
      }
    }

  ifs.close();
  
}


void ConvertTrackingData(TrackingDataList& coordinates, MatrixList& matrices)
{
  
  TrackingDataList::iterator iter;

  for (iter = coordinates.begin(); iter != coordinates.end(); iter ++)
    {
    MatrixData matrixData;
    matrixData.ts = iter->ts;
    
    igtl::Matrix4x4& matrix1 = matrixData.matrixCoil1;
    igtl::Matrix4x4& matrix2 = matrixData.matrixCoil2;
    igtl::Matrix4x4& matrix3 = matrixData.matrixCoil3;
    igtl::Matrix4x4& matrix4 = matrixData.matrixCoil4;
    
    igtl::IdentityMatrix(matrix1);
    igtl::IdentityMatrix(matrix2);
    igtl::IdentityMatrix(matrix3);
    igtl::IdentityMatrix(matrix4);
    
    matrixData.ts = iter->ts;
    matrix1[0][3] = iter->coil1X;
    matrix1[1][3] = iter->coil1Y;
    matrix1[2][3] = iter->coil1Z;
    matrix2[0][3] = iter->coil2X;
    matrix2[1][3] = iter->coil2Y;
    matrix2[2][3] = iter->coil2Z;
    matrix3[0][3] = iter->coil3X;
    matrix3[1][3] = iter->coil3Y;
    matrix3[2][3] = iter->coil3Z;
    matrix4[0][3] = iter->coil4X;
    matrix4[1][3] = iter->coil4Y;
    matrix4[2][3] = iter->coil4Z;

    matrices.push_back(matrixData);
    }
}

int SendTrackingData(igtl::ClientSocket::Pointer& socket, igtl::TrackingDataMessage::Pointer& trackingMsg, MatrixData& mat)
{
  igtl::Matrix4x4 matrix;
  igtl::TrackingDataElement::Pointer ptr;

  std::cout << "===== Time: " << mat.ts << " =====" << std::endl;

  // Coil 1
  trackingMsg->GetTrackingDataElement(0, ptr);
  ptr->SetMatrix(mat.matrixCoil1);
  igtl::PrintMatrix(mat.matrixCoil1);
  
  // Coil 2
  trackingMsg->GetTrackingDataElement(1, ptr);
  ptr->SetMatrix(mat.matrixCoil2);
  igtl::PrintMatrix(mat.matrixCoil2);
  
  // Coil 3
  trackingMsg->GetTrackingDataElement(2, ptr);
  ptr->SetMatrix(mat.matrixCoil3);
  igtl::PrintMatrix(mat.matrixCoil3);

  // Coil 4
  trackingMsg->GetTrackingDataElement(3, ptr);
  ptr->SetMatrix(mat.matrixCoil4);
  igtl::PrintMatrix(mat.matrixCoil4);

  trackingMsg->Pack();
  socket->Send(trackingMsg->GetPackPointer(), trackingMsg->GetPackSize());
  
  return 0;
}
  
int SendTrackingData(igtl::ClientSocket::Pointer& socket, igtl::TrackingDataMessage::Pointer& trackingMsg)
{

  static float phi0   = 0.0;
  static float theta0 = 0.0;
  static float phi1   = 0.0;
  static float theta1 = 0.0;
  static float phi2   = 0.0;
  static float theta2 = 0.0;

  igtl::Matrix4x4 matrix;
  igtl::TrackingDataElement::Pointer ptr;

  // Coil 0
  trackingMsg->GetTrackingDataElement(0, ptr);
  GetRandomTestMatrix(matrix, phi0, theta0);
  ptr->SetMatrix(matrix);
  
  // Coil 1
  trackingMsg->GetTrackingDataElement(1, ptr);
  GetRandomTestMatrix(matrix, phi1, theta1);
  ptr->SetMatrix(matrix);
  
  // Coil 2
  trackingMsg->GetTrackingDataElement(2, ptr);
  GetRandomTestMatrix(matrix, phi2, theta2);
  ptr->SetMatrix(matrix);

  trackingMsg->Pack();
  socket->Send(trackingMsg->GetPackPointer(), trackingMsg->GetPackSize());
  
  phi0 += 0.1;
  phi1 += 0.2;
  phi2 += 0.3;
  theta0 += 0.2;
  theta1 += 0.1;
  theta2 += 0.05;

  return 0;
}



//------------------------------------------------------------
// Function to generate random matrix.
void GetRandomTestMatrix(igtl::Matrix4x4& matrix, float phi, float theta)
{
  float position[3];
  float orientation[4];

  // random position
  position[0] = 50.0 * cos(phi);
  position[1] = 50.0 * sin(phi);
  position[2] = 50.0 * cos(phi);
  phi = phi + 0.2;

  // random orientation
  orientation[0]=0.0;
  orientation[1]=0.6666666666*cos(theta);
  orientation[2]=0.577350269189626;
  orientation[3]=0.6666666666*sin(theta);
  theta = theta + 0.1;

  //igtl::Matrix4x4 matrix;
  igtl::QuaternionToMatrix(orientation, matrix);

  matrix[0][3] = position[0];
  matrix[1][3] = position[1];
  matrix[2][3] = position[2];
  
  //igtl::PrintMatrix(matrix);
}






