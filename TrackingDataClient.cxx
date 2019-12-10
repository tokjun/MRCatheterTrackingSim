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
#include "igtlTransformMessage.h"
#include "igtlClientSocket.h"
#include "igtlMath.h"

typedef struct {
  double ts;
  double V0X;  //  tip1X 
  double V0Y;  //  tip1Y 
  double V0Z;  //  tip1Z 
  double V1X;  //  tip1NX
  double V1Y;  //  tip1NY
  double V1Z;  //  tip1NZ
  double V2X;  //  coil1X
  double V2Y;  //  coil1Y
  double V2Z;  //  coil1Z
  double V3X;  //  coil2X
  double V3Y;  //  coil2Y
  double V3Z;  //  coil2Z
  double V4X;  //  tip2X 
  double V4Y;  //  tip2Y 
  double V4Z;  //  tip2Z 
  double V5X;  //  tip2NX
  double V5Y;  //  tip2NY
  double V5Z;  //  tip2NZ
  double V6X;  //  coil3X
  double V6Y;  //  coil3Y
  double V6Z;  //  coil3Z
  double V7X;  //  coil4X
  double V7Y;  //  coil4Y
  double V7Z;  //  coil4Z
} TrackingData;


typedef struct {
  double ts;
  igtl::Matrix4x4 matrix;
} MatrixFrame;


typedef struct {
  double ts;
  igtl::Matrix4x4 matrixCoil1;
  igtl::Matrix4x4 matrixCoil2;
  igtl::Matrix4x4 matrixCoil3;
  igtl::Matrix4x4 matrixCoil4;
} CatheterMatricesFrame;

enum {
  TRANSFORM,
  TRACKING,
};

typedef std::vector<TrackingData> TrackingDataList;
typedef std::vector<MatrixFrame> MatrixFrameList;
typedef std::vector<CatheterMatricesFrame> CatheterMatricesFrameList;


void  ReadFile(std::string filename, TrackingDataList& coordinates);
void  ConvertTransformData(TrackingDataList& coordinates, MatrixFrameList& frameList);
void  ConvertTrackingData(TrackingDataList& coordinates, CatheterMatricesFrameList& frameList);
int   SendTransformData(igtl::ClientSocket::Pointer& socket, igtl::TransformMessage::Pointer& transformMsg, MatrixFrame& mat);
int   SendTrackingData(igtl::ClientSocket::Pointer& socket, igtl::TrackingDataMessage::Pointer& trackingMsg, CatheterMatricesFrame& mat);

int main(int argc, char* argv[])
{
  //------------------------------------------------------------
  // Parse Arguments

  if (argc != 6) // check number of arguments
    {
    // If not correct, print usage
    std::cerr << "Usage: " << argv[0] << " <hostname> <port> <IGTL type> <fps> <file>"    << std::endl;
    std::cerr << "    <hostname>  : IP or host name"                    << std::endl;
    std::cerr << "    <port>      : Port # (18944 in Slicer default)"   << std::endl;
    std::cerr << "    <IGTL type> : M = Transform; T = Tracking data"   << std::endl;
    std::cerr << "    <fps>       : Frequency (fps) to send coordinate" << std::endl;
    std::cerr << "    <file>      : Tracking file" << std::endl;
    exit(0);
    }

  char*  hostname = argv[1];
  int    port     = atoi(argv[2]);
  int    type;
  if (argv[3][0] == 'M')
    {
    type = TRANSFORM;
    }
  else
    {
    type = TRACKING;
    }
  double fps      = atof(argv[4]);
  int    interval = (int) (1000.0 / fps);
  std::string filename = argv[5];


  //------------------------------------------------------------
  // Load Tracking Data
  
  TrackingDataList coordinates;
  MatrixFrameList mFrameList;
  CatheterMatricesFrameList cmFrameList;

  ReadFile(filename, coordinates);
  
  if (type == TRANSFORM)
    {
    ConvertTransformData(coordinates, mFrameList);
    }
  else
    {
    ConvertTrackingData(coordinates, cmFrameList);
    }
  
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

  if (type == TRANSFORM)
    {
    //------------------------------------------------------------
    // Allocate Transform Message Class
    //
    
    igtl::TransformMessage::Pointer transMsg;
    transMsg = igtl::TransformMessage::New();
    
    //------------------------------------------------------------
    // Loop
    MatrixFrameList::iterator iter;
    for (iter = mFrameList.begin(); iter != mFrameList.end(); iter ++)
      {
      transMsg->SetDeviceName("Tracker");
      SendTransformData(socket, transMsg, *iter);
      igtl::Sleep(interval);
      }

    }
  else // type = TRACKING
    {
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
    CatheterMatricesFrameList::iterator iter;
    for (iter = cmFrameList.begin(); iter != cmFrameList.end(); iter ++)
      {
      trackingMsg->SetDeviceName("Tracker");
      SendTrackingData(socket, trackingMsg, *iter);
      igtl::Sleep(interval);
      }
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
        >> pt.V0X >> pt.V0Y >> pt.V0Z
        >> pt.V1X >> pt.V1Y >> pt.V1Z
        >> pt.V2X >> pt.V2Y >> pt.V2Z
        >> pt.V3X >> pt.V3Y >> pt.V3Z
        >> pt.V4X >> pt.V4Y >> pt.V4Z
        >> pt.V5X >> pt.V5Y >> pt.V5Z
        >> pt.V6X >> pt.V6Y >> pt.V6Z
        >> pt.V7X >> pt.V7Y >> pt.V7Z)
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


void  ConvertTransformData(TrackingDataList& coordinates, MatrixFrameList& frameList)
{

  frameList.clear();
  
  // Assuming that the first and second vectors (V0 and V1) represent the tip position
  // and the catheter orientation respectively.

  TrackingDataList::iterator iter;

  for (iter = coordinates.begin(); iter != coordinates.end(); iter ++)
    {
    MatrixFrame frame;
    igtl::Matrix4x4& matrix = frame.matrix;
  
    float t[3];
    float s[3];
    float n[3];
    float nlen;
    
    frame.ts = iter->ts;
    
    n[0] = iter->V1X;
    n[1] = iter->V1Y;
    n[2] = iter->V1Z;
    
    nlen = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (nlen > 0)
      {
      n[0] /= nlen;
      n[1] /= nlen;
      n[2] /= nlen;
      }
    else
      {
      n[0] = 0.0;
      n[1] = 0.0;
      n[2] = 1.0;
      }
    
    // Check if <n> is not parallel to <s>=(0.0, 1.0, 0.0)
    if (n[1] < 1.0)
      {
      s[0] = 0.0;
      s[1] = 1.0;
      s[2] = 0.0;
      igtl::Cross(t, s, n);
      igtl::Cross(s, n, t);
      }
    else
      {
      t[0] = 1.0;
      t[1] = 0.0;
      t[2] = 0.0;
      igtl::Cross(s, n, t);
      igtl::Cross(t, s, n);
      }
    
    matrix[0][0] = t[0];
    matrix[1][0] = t[1];
    matrix[2][0] = t[2];
    matrix[0][1] = s[0];
    matrix[1][1] = s[1];
    matrix[2][1] = s[2];
    matrix[0][2] = n[0];
    matrix[1][2] = n[1];
    matrix[2][2] = n[2];
    matrix[0][3] = iter->V0X;
    matrix[1][3] = iter->V0Y;
    matrix[2][3] = iter->V0Z;
    
    igtl::PrintMatrix(matrix);
    frameList.push_back(frame);
    }
    
}


void  ConvertTrackingData(TrackingDataList& coordinates, CatheterMatricesFrameList& frameList)
{

  frameList.clear();
  
  TrackingDataList::iterator iter;

  for (iter = coordinates.begin(); iter != coordinates.end(); iter ++)
    {
    CatheterMatricesFrame frame;
    frame.ts = iter->ts;
    
    igtl::Matrix4x4& matrix1 = frame.matrixCoil1;
    igtl::Matrix4x4& matrix2 = frame.matrixCoil2;
    igtl::Matrix4x4& matrix3 = frame.matrixCoil3;
    igtl::Matrix4x4& matrix4 = frame.matrixCoil4;
    
    igtl::IdentityMatrix(matrix1);
    igtl::IdentityMatrix(matrix2);
    igtl::IdentityMatrix(matrix3);
    igtl::IdentityMatrix(matrix4);

    frame.ts = iter->ts;
    matrix1[0][3] = iter->V2X;
    matrix1[1][3] = iter->V2Y;
    matrix1[2][3] = iter->V2Z;
    matrix2[0][3] = iter->V3X;
    matrix2[1][3] = iter->V3Y;
    matrix2[2][3] = iter->V3Z;
    matrix3[0][3] = iter->V6X;
    matrix3[1][3] = iter->V6Y;
    matrix3[2][3] = iter->V6Z;
    matrix4[0][3] = iter->V7X;
    matrix4[1][3] = iter->V7Y;
    matrix4[2][3] = iter->V7Z;

    frameList.push_back(frame);
    }
}


int SendTransformData(igtl::ClientSocket::Pointer& socket, igtl::TransformMessage::Pointer& transformMsg, MatrixFrame& mat)
{
  igtl::Matrix4x4 matrix;

  std::cout << "===== Time: " << mat.ts << " =====" << std::endl;
  igtl::PrintMatrix(mat.matrix);

  transformMsg->SetMatrix(mat.matrix);
  transformMsg->Pack();
  socket->Send(transformMsg->GetPackPointer(), transformMsg->GetPackSize());
  
  return 0;
}


int SendTrackingData(igtl::ClientSocket::Pointer& socket, igtl::TrackingDataMessage::Pointer& trackingMsg, CatheterMatricesFrame& mat)
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




