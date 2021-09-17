/*=========================================================================

  Program:   MR Tracking Simulator
  Language:  C++

  Copyright (c) Junichi Tokuda

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
#include <vector>
#include <iomanip>

#include "igtlOSUtil.h"
#include "igtlTrackingDataMessage.h"
#include "igtlTransformMessage.h"
#include "igtlStringMessage.h"
#include "igtlClientSocket.h"
#include "igtlMath.h"


#define MAX_CHANNELS 8

enum FrameType { TYPE_TRACKING, TYPE_TRANSFORM, TYPE_STRING };

// Vector (x, y, z)
typedef struct {
  double x;
  double y;
  double z;
} Vector;

// 4x4 matrix (Upper 3 rows only assuming that the bottom row is (0, 0, 0, 1))
typedef struct {
  float m00;
  float m01;
  float m02;
  float m03;
  float m10;
  float m11;
  float m12;
  float m13;
  float m20;
  float m21;  
  float m22;  
  float m23;
  //float m30;
  //float m31;
  //float m32;
  //float m33;  
} Matrix4x4;

typedef struct {
  double ts;
  int type;   // FrameType
  std::string name;
  std::vector<Vector> vectors;
  std::string text;
} RawData;

typedef struct {
  double ts;
  int type;  // FrameType
  std::string name;
  std::vector<Matrix4x4> matrices;
  std::string text;
} Frame;

enum {
  TRANSFORM,
  TRACKING,
};

enum {
  POSITION,
  NORMAL,
};

typedef std::vector<RawData> RawDataList;
typedef std::vector<Frame> FrameList;

void  ReadFile(std::string filename, RawDataList& coordinates, char delim, bool fAscending);
void  ConvertTransformData(RawDataList& coordinates, FrameList& frameList, bool fNormal);
void  ConvertTrackingData(RawDataList& coordinates, FrameList& frameList);
int   SendTransformData(igtl::ClientSocket::Pointer& socket, igtl::TransformMessage::Pointer& transformMsg, Frame& mat);
int   SendTrackingData(igtl::ClientSocket::Pointer& socket, igtl::TrackingDataMessage::Pointer& trackingMsg, Frame& mat, std::vector<bool> mask);
int   SendStringData(igtl::ClientSocket::Pointer& socket, igtl::StringMessage::Pointer& stringMsg, Frame& mat);


void SplitString(char* input, std::vector<std::string>& output, char delim, int len)
{
  // Split the current line by delimiter
  char* p = input;
  char* pend = p + len;
  
  output.clear();
  while(p < pend && *p != '\n')
    {
    char buf[1024];
    char* pbuf = buf;
    while (*p != delim && p < pend)
      {
      *pbuf = *p;
      p ++;
      pbuf ++;
      }
    *pbuf = '\0';
    
    // Ignore following delimiters
    while (*p == delim && p < pend)
      {
      p++;
      }
    std::string col = buf;
    output.push_back(col);
    }
}


void PrintUsage(const char* progName)
{
  std::cerr << "Usage: " << progName << " [-h <hostname>] [-p <port>] [-t <IGTL type>] [-f <fps>] [-m <mask>] [-n] [-o <order>] [-D] [-d <dev name>] <file>"    << std::endl;
  std::cerr << "    -h <hostname>  : IP or host name (default: \"localhost\")" << std::endl;
  std::cerr << "    -p <port>      : Port # (default: \"18944\")"   << std::endl;
  std::cerr << "    -t <IGTL type> : Output messagr type. 'T' = Tracking data (default); 'M' = Transform; when 'M' is specified, only the first two sets of vectors are used." << std::endl;
  std::cerr << "    -f <fps>       : Frequency (fps) to send coordinate (default: 5); if 't' is specified, the interval will be calculated from time stamp. " << std::endl;
  std::cerr << "    -m <mask>      : Channel mask (Maximum " << MAX_CHANNELS << " channels); ex. '11111111', '10101010' (default: \"11111111\")" << std::endl;
  std::cerr << "    -n             : Specify if normal vector is used." << std::endl;
  std::cerr << "    -o <order>     : Order of channels ('a' for ascending / 'd' for desending; default is 'a')" << std::endl;
  std::cerr << "    -d <dev name>  : Device name (default \"Tracking\")" << std::endl;
  std::cerr << "    <file>      : Tracking file" << std::endl;
}

int main(int argc, char* argv[])
{
  //------------------------------------------------------------
  // Parse Arguments

  // Variables and default values
  std::string hostname       = "localhost";
  int    port                = 18944;
  int    type                = TRACKING;
  bool   fNormal             = false;
  bool   fAscending          = true;
  double fps                 = 5.0;
  bool   fUseTimeStamp       = false;
  std::string devName        = "Tracking";
  std::vector<bool> chmask;
  std::string filename = "";

  chmask.clear();
  
  for (int i = 1; i < argc; i ++)
    {
    if (strncmp(argv[i], "-h", 2) == 0 && i < argc-1)
      {
      i ++;
      hostname = argv[i];
      }
    else if (strncmp(argv[i], "-p", 2) == 0 && i < argc-1)
      {
      i ++;
      port = atoi(argv[i]);
      }
    else if (strncmp(argv[i], "-t", 2) == 0 && i < argc-1)
      {
      i ++;
      if (argv[i][0] == 'M')
        {
        type = TRANSFORM;
        }
      else
        {
        type = TRACKING;
        }
      }
    else if (strncmp(argv[i], "-n", 2) == 0 && i < argc-1)
      {
      fNormal = true;
      }
    else if (strncmp(argv[i], "-o", 2) == 0 && i < argc-1)
      {
      i ++;
      if (argv[i][0] == 'd')
        {
        fAscending = false;
        }
      else if (argv[i][0] == 'a')
        {
        fAscending = true;
        }
      else
        {
        std::cerr << "Wrong option for -o: " << argv[i] << std::endl;
        exit(0);
        }
      }
    else if (strncmp(argv[i], "-f", 2) == 0 && i < argc-1)
      {
      i ++;
      if (strncmp(argv[i], "t", 1) == 0)
        {
        fUseTimeStamp = true;
        }
      else
        {
        fps = atof(argv[i]);
        }
      }
    else if (strncmp(argv[i], "-d", 2) == 0 && i < argc-1)
      {
      i ++;      
      devName  = argv[i];
      }
    else if (strncmp(argv[i], "-m", 2) == 0 && i < argc-1)
      {
      i ++;
      if (i < argc - 1 && strlen(argv[i]) <= MAX_CHANNELS)
        {
        for (int j = 0; j < strlen(argv[i]); j ++)
          {
            if (argv[i][j] == '1')
            {
            chmask.push_back(true);
            }
          else
            {
            chmask.push_back(false);
            }
          }
        }
      else
        {
        std::cerr << "The mask size exceeded the limit (" << MAX_CHANNELS << ")" << std::endl;
        exit(0);
        }
      }
    else 
      {
      if (strlen(argv[i]) == 2 && argv[i][0] == '-')
        {
        // invalid parameter
        std::cerr << "Invalid parameter: " << argv[i] << std::endl << std::endl;
        PrintUsage(argv[0]);
        exit(0);
        }
      else // Suppose file name parameter
        {
        if (filename.compare("") != 0)
          {
          std::cerr << "Invalid parameter: " << argv[i] << std::endl << std::endl;
          PrintUsage(argv[0]);
          exit(0);
          }
        else
          {
          filename = argv[i];
          }
        }
      }
    }
    
  int    interval = (int) (1000.0 / fps);

  if (filename.compare("") == 0)
    {
    PrintUsage(argv[0]);
    exit(0);
    }

  if (chmask.size() == 0)
    {
    for (int i = 0; i < MAX_CHANNELS; i ++)
      {
      chmask.push_back(true);
      }
    }

  int nCh = 0;
  std::vector<bool>::iterator iter;
  for (iter = chmask.begin(); iter != chmask.end(); iter ++)
    {
    if (*iter) nCh ++;
    }

  std::cout << "Number of channels = " << nCh << "/" << chmask.size() << std::endl;

  //------------------------------------------------------------
  // Load Tracking Data
  
  RawDataList coordinates;
  FrameList mFrameList;
  //CatheterMatricesFrameList cmFrameList;

  ReadFile(filename, coordinates, '\t', fAscending);
  
  if (type == TRANSFORM)
    {
    ConvertTransformData(coordinates, mFrameList, fNormal);
    }
  else
    {
    ConvertTrackingData(coordinates, mFrameList);
    }
  
  //------------------------------------------------------------
  // Establish Connection

  igtl::ClientSocket::Pointer socket;
  socket = igtl::ClientSocket::New();
  int r = socket->ConnectToServer(hostname.c_str(), port);

  if (r != 0)
    {
    std::cerr << "Cannot connect to the server." << std::endl;
    exit(0);
    }

  //------------------------------------------------------------
  // Allocate Transform Message Class
  //
  
  igtl::TransformMessage::Pointer transMsg;
  transMsg = igtl::TransformMessage::New();
    
  //------------------------------------------------------------
  // Allocate TrackingData Message Class
  //
  // NOTE: TrackingDataElement class instances are allocated
  //       before the loop starts to avoid reallocation
  //       in each image transfer.
  
  igtl::TrackingDataMessage::Pointer trackingMsg;
  trackingMsg = igtl::TrackingDataMessage::New();
  
  igtl::TrackingDataElement::Pointer * trackElement;
  trackElement = new igtl::TrackingDataElement::Pointer[nCh];
  
  for (int i = 0; i < nCh; i ++)
    {
    std::stringstream ss;
    ss << devName.c_str() << std::setfill('0') << std::setw(2) << i;
    trackElement[i] = igtl::TrackingDataElement::New();
    trackElement[i]->SetName(ss.str().c_str());
    trackElement[i]->SetType(igtl::TrackingDataElement::TYPE_3D);
    trackingMsg->AddTrackingDataElement(trackElement[i]);
    }

  //------------------------------------------------------------
  // Allocate Transform Message Class

  igtl::StringMessage::Pointer stringMsg;
  stringMsg = igtl::StringMessage::New();

  //------------------------------------------------------------
  // Loop
  FrameList::iterator fliter;
  double prevTimeStamp = -1;
  
  igtl::TimeStamp::Pointer cts;
  cts = igtl::TimeStamp::New();
  double stt;
  double elapsed = 0.0;
  
  for (fliter = mFrameList.begin(); fliter != mFrameList.end(); fliter ++)
    {
      
    // Wait
    if (fUseTimeStamp)
      {
      int _interval; 
      if (prevTimeStamp < 0)
        {
        _interval = 0;
        }
      else
        {
        _interval = (int) ((fliter->ts - prevTimeStamp - elapsed) * 1000.0);
        if (_interval < 0)
          {
          _interval = 0;
          }
        }
      std::cerr << fliter->type << " Time: " << elapsed << " / " << _interval << std::endl;
      igtl::Sleep(_interval);
      prevTimeStamp = fliter->ts;
      }
    else
      {
      igtl::Sleep(interval);
      }
    
    cts->GetTime();
    stt = cts->GetTimeStamp();
    
    // Set up meesage
    if (fliter->type == TYPE_TRANSFORM)
      {
      transMsg->SetDeviceName(fliter->name.c_str());
      SendTransformData(socket, transMsg, *fliter);
      }
    else if (fliter->type == TYPE_TRACKING)
      {
      trackingMsg->SetDeviceName(fliter->name.c_str());
      SendTrackingData(socket, trackingMsg, *fliter, chmask);
      }
    else // String
      {
      stringMsg->SetDeviceName(fliter->name.c_str());
      SendStringData(socket, stringMsg, *fliter);
      }
    cts->GetTime();
    elapsed = cts->GetTimeStamp() - stt;
    
    }
  delete[] trackElement;
      
}


void  ReadFile(std::string filename, RawDataList& coordinates, char delim, bool fAscending)
{
  
  std::ifstream ifs;

  ifs.open(filename);
  if (!ifs.is_open())
    {
    std::cerr << "Error opening file: " << filename << std::endl;
    exit(0);
    }

  char line[1024];

  // Parse lines
  while (ifs.getline(line, 1023))
    {
    if (ifs.eof())
      {
      break;
      }
    
    std::vector<std::string> cols;
    int len = strnlen(line, 1023);
    SplitString(line, cols, delim, len);
    
    RawData pt;
    pt.vectors.clear();
    
    pt.name = cols[1].c_str();
    pt.ts = atof(cols[2].c_str());
    //std::cerr << cols[0] << " - "  << pt.name << std::endl;
    
    // Type
    std::string typeStr = cols[0];
    int offset=3;
    if (typeStr.compare("TDATA") == 0)
      {
      int nvec; // Number of data columns in the line
      nvec = (cols.size() - offset) / 3;
      pt.type = TYPE_TRACKING;
      for (int i = 0; i < nvec; i ++)
        {
        Vector v;
        v.x = atof(cols[offset+i*3+0].c_str());
        v.y = atof(cols[offset+i*3+1].c_str());
        v.z = atof(cols[offset+i*3+2].c_str());
        if (fAscending)
          {
          pt.vectors.push_back(v);
          }
        else
          {
          pt.vectors.insert(pt.vectors.begin(), v);
          }
        }
      }
    else // typeStr == "STRING"
      {
      int nlines; // Number of columns (=lines) for string data
      nlines = cols.size() - 3;
      pt.type = TYPE_STRING;
      pt.text = "";
      for (int i = 0; i < nlines; i ++)
        {
        pt.text += cols[offset+i] + "\n";
        }
      }
    // std::cerr << "Read: " << coordinates.size() << " coordinates from file." << std::endl;
    std::cerr << "ReadFile  " << pt.name << " " << pt.type << std::endl;
    coordinates.push_back(pt);
    }
  
  ifs.close();
}


void  ConvertTransformData(RawDataList& coordinates, FrameList& frameList, bool fNormal)
{

  frameList.clear();

  // If fNormal=true, we assume that the first and second vectors (V0 and V1) represent the tip position
  // and the catheter orientation respectively

  RawDataList::iterator iter;

  for (iter = coordinates.begin(); iter != coordinates.end(); iter ++)
    {
    Frame frame;
    frame.ts = iter->ts;
    frame.name = iter->name;
    
    //std::cerr << "ConvertTransformData  " << frame.name << " " << frame.type << std::endl;
    
    if (iter->type == TYPE_TRACKING)
      {
      frame.type = TYPE_TRANSFORM;
      frame.matrices.resize(1);
      Matrix4x4& matrix = frame.matrices[0];
      
      float t[3];
      float s[3];
      float n[3];
      float nlen;
      
      // Make sure to have more than two vectors in the tracking data
      
      if (fNormal)
        {
        n[0] = iter->vectors[1].x;
        n[1] = iter->vectors[1].y;
        n[2] = iter->vectors[1].z;
        }
      else
        {
        n[0] = iter->vectors[0].x - iter->vectors[1].x;
        n[1] = iter->vectors[0].y - iter->vectors[1].y;
        n[2] = iter->vectors[0].z - iter->vectors[1].z;
        }
        
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
      
      matrix.m00 = t[0];
      matrix.m10 = t[1];
      matrix.m20 = t[2];
      matrix.m01 = s[0];
      matrix.m11 = s[1];
      matrix.m21 = s[2];
      matrix.m02 = n[0];
      matrix.m12 = n[1];
      matrix.m22 = n[2];
      matrix.m03 = iter->vectors[0].x;
      matrix.m13 = iter->vectors[0].y;
      matrix.m23 = iter->vectors[0].z;
      
      //igtl::PrintMatrix(matrix);
      }
    else // iter->type == TYPE_STRING
      {
      frame.type = TYPE_STRING;
      frame.text = iter->text;
      }
    frameList.push_back(frame);
    }
}


void  ConvertTrackingData(RawDataList& coordinates, FrameList& frameList)
{

  frameList.clear();
  
  RawDataList::iterator iter;

  for (iter = coordinates.begin(); iter != coordinates.end(); iter ++)
    {
    Frame frame;
    frame.ts = iter->ts;
    frame.name = iter->name;
    //std::cerr << "ConverTrackingData()  " << frame.name << " " << frame.type << std::endl;
    
    if (iter->type == TYPE_TRACKING)
      {
      frame.type = TYPE_TRACKING;
      //Current size
      int len = iter->vectors.size();
      frame.type = TYPE_TRACKING;
      frame.matrices.resize(len);
      
      for (int i = 0; i < len; i ++)
        {
        Matrix4x4& mat = frame.matrices[i];
        mat.m00 = 1.0;
        mat.m10 = 0.0;
        mat.m20 = 0.0;
        mat.m01 = 0.0;
        mat.m11 = 1.0;
        mat.m21 = 0.0;
        mat.m02 = 0.0;
        mat.m12 = 0.0;
        mat.m22 = 1.0;
        mat.m03 = iter->vectors[i].x;
        mat.m13 = iter->vectors[i].y;
        mat.m23 = iter->vectors[i].z;
        }
    
      }
    else // iter->type == TYPE_STRING
      {
      frame.type = TYPE_STRING;
      frame.text = iter->text;
      }
    frameList.push_back(frame);
    }

}


int SendTransformData(igtl::ClientSocket::Pointer& socket, igtl::TransformMessage::Pointer& transformMsg, Frame& mat)
{
  if (mat.matrices.size() < 1)
    {
    return 0;
    }
  
  std::cout << "===== Time: " << mat.ts << " / " << mat.name << " =====" << std::endl;
  igtl::Matrix4x4 matrix;
  Matrix4x4& m = mat.matrices[0];
  matrix[0][0] = m.m00;
  matrix[0][1] = m.m01;
  matrix[0][2] = m.m02;
  matrix[0][3] = m.m03;
  matrix[1][0] = m.m10;
  matrix[1][1] = m.m11;
  matrix[1][2] = m.m12;
  matrix[1][3] = m.m13;
  matrix[2][0] = m.m20;
  matrix[2][1] = m.m21;
  matrix[2][2] = m.m22;
  matrix[2][3] = m.m23;
  igtl::PrintMatrix(matrix);

  transformMsg->SetMatrix(matrix);
  transformMsg->Pack();
  socket->Send(transformMsg->GetPackPointer(), transformMsg->GetPackSize());
  
  return 1;
}


int SendTrackingData(igtl::ClientSocket::Pointer& socket, igtl::TrackingDataMessage::Pointer& trackingMsg, Frame& mat, std::vector<bool> mask)
{
  igtl::Matrix4x4 matrix;
  igtl::TrackingDataElement::Pointer ptr;

  igtl::IdentityMatrix(matrix);
  
  int nCh = 0;
  std::vector<bool>::iterator iter;
  for (iter = mask.begin(); iter != mask.end(); iter ++)
    {
    if (*iter) nCh ++;
    }

  if (mat.matrices.size() < mask.size())
    {
    std::cerr << "Mask size is larger than the number of tracking channels." << std::endl;
    return 0;
    }
  
  std::cout << "===== Time: " << mat.ts << " / " << mat.name << " =====" << std::endl;

  std::vector<bool>::iterator mskIter;
  std::vector<Matrix4x4>::iterator matIter;

  int ch = 0;
  matIter = mat.matrices.begin();
  
  for (mskIter = mask.begin(); mskIter != mask.end(); mskIter ++)
    {
    if (*mskIter) // Mask is on
      {
      matrix[0][0] = matIter->m00;
      matrix[0][1] = matIter->m01;
      matrix[0][2] = matIter->m02;
      matrix[0][3] = matIter->m03;
      matrix[1][0] = matIter->m10;
      matrix[1][1] = matIter->m11;
      matrix[1][2] = matIter->m12;
      matrix[1][3] = matIter->m13;
      matrix[2][0] = matIter->m20;
      matrix[2][1] = matIter->m21;
      matrix[2][2] = matIter->m22;
      matrix[2][3] = matIter->m23;
                               
      trackingMsg->GetTrackingDataElement(ch, ptr);
      ptr->SetMatrix(matrix);
      igtl::PrintMatrix(matrix);
      ch ++;
      }
    matIter ++;
    }
    
  trackingMsg->Pack();
  socket->Send(trackingMsg->GetPackPointer(), trackingMsg->GetPackSize());
  
  return ch;
}

int   SendStringData(igtl::ClientSocket::Pointer& socket, igtl::StringMessage::Pointer& stringMsg, Frame& mat)
{
  std::cerr << "Sending String" << std::endl;
  stringMsg->SetString(mat.text);
  stringMsg->Pack();
  socket->Send(stringMsg->GetPackPointer(), stringMsg->GetPackSize());

  return 1;
}


