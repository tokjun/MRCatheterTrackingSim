/*=========================================================================

  Program:   OpenIGTLink -- Example for Tracker Client Program
  Language:  C++

  Copyright (c) Insight Software Consortium. All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <iostream>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <iomanip>

#include "igtlOSUtil.h"
#include "igtlTrackingDataMessage.h"
#include "igtlClientSocket.h"

int   SendTrackingData(igtl::ClientSocket::Pointer& socket, igtl::TrackingDataMessage::Pointer& trackingMsg);
void  GetRandomTestMatrix(igtl::Matrix4x4& matrix, float phi, float theta);

int main(int argc, char* argv[])
{
  //------------------------------------------------------------
  // Parse Arguments

  if (argc != 4) // check number of arguments
    {
    // If not correct, print usage
    std::cerr << "Usage: " << argv[0] << " <hostname> <port> <fps>"    << std::endl;
    std::cerr << "    <hostname> : IP or host name"                    << std::endl;
    std::cerr << "    <port>     : Port # (18944 in Slicer default)"   << std::endl;
    std::cerr << "    <fps>      : Frequency (fps) to send coordinate" << std::endl;
    exit(0);
    }

  char*  hostname = argv[1];
  int    port     = atoi(argv[2]);
  double fps      = atof(argv[3]);
  int    interval = (int) (1000.0 / fps);

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

  igtl::TrackingDataElement::Pointer trackElement[3];
  
  for (int coil = 0; coil < 3; coil ++)
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
  while (1)
    {
    trackingMsg->SetDeviceName("Tracker");
    SendTrackingData(socket, trackingMsg);
    igtl::Sleep(interval);
    }
      
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

  // Channel 0
  trackingMsg->GetTrackingDataElement(0, ptr);
  GetRandomTestMatrix(matrix, phi0, theta0);
  ptr->SetMatrix(matrix);
  
  // Channel 1
  trackingMsg->GetTrackingDataElement(1, ptr);
  GetRandomTestMatrix(matrix, phi1, theta1);
  ptr->SetMatrix(matrix);
  
  // Channel 2
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






