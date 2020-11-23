#include "trojanmap.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <algorithm>
#include <fstream>
#include <locale>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <stack>

#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/videoio.hpp"

//-----------------------------------------------------
// TODO (Students): You do not and should not change the following functions:
//-----------------------------------------------------

/**
 * PrintMenu: Create the menu
 * 
 */
void TrojanMap::PrintMenu() {

  std::string menu =
      "**************************************************************\n"
      "* Select the function you want to execute.                    \n"
      "* 1. Autocomplete                                             \n"
      "* 2. Find the position                                        \n"
      "* 3. CalculateShortestPath                                    \n"
      "* 4. Travelling salesman problem                              \n"
      "* 5. Exit                                                     \n"
      "**************************************************************\n";
  std::cout << menu << std::endl;
  std::string input;
  getline(std::cin, input);
  char number = input[0];
  switch (number)
  {
  case '1':
  {
    menu =
        "**************************************************************\n"
        "* 1. Autocomplete                                             \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "Please input a partial location:";
    std::cout << menu;
    getline(std::cin, input);
    auto results = Autocomplete(input);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    if (results.size() != 0) {
      for (auto x : results) std::cout << x << std::endl;
    } else {
      std::cout << "No matched locations." << std::endl;
    }
    menu = "**************************************************************\n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }
  case '2':
  {
    menu =
        "**************************************************************\n"
        "* 2. Find the position                                        \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "Please input a location:";
    std::cout << menu;
    getline(std::cin, input);
    auto results = GetPosition(input);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    if (results.first != -1) {
      std::cout << "Latitude: " << results.first
                << " Longitude: " << results.second << std::endl;
      PlotPoint(results.first, results.second);
    } else {
      std::cout << "No matched locations." << std::endl;
    }
    menu = "**************************************************************\n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }
  case '3':
  {
    menu =
        "**************************************************************\n"
        "* 3. CalculateShortestPath                                            "
        "      \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "Please input the start location:";
    std::cout << menu;
    std::string input1;
    getline(std::cin, input1);
    menu = "Please input the destination:";
    std::cout << menu;
    std::string input2;
    getline(std::cin, input2);
    auto results = CalculateShortestPath(input1, input2);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    if (results.size() != 0) {
      for (auto x : results) std::cout << x << std::endl;
      PlotPath(results);
    } else {
      std::cout << "No route from the start point to the destination."
                << std::endl;
    }
    menu = "**************************************************************\n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }
  case '4':
  {
    menu =
        "**************************************************************\n"
        "* 4. Travelling salesman problem                              \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "In this task, we will select N random points on the map and you need to find the path to travel these points and back to the start point.";
    std::cout << menu << std::endl << std::endl;
    menu = "Please input the number of the places:";
    std::cout << menu;
    getline(std::cin, input);
    int num = std::stoi(input);
    std::vector<std::string> keys;
    for (auto x : data) {
      keys.push_back(x.first);
    }
    std::vector<std::string> locations;
    srand(time(NULL));
    for (int i = 0; i < num; i++)
      locations.push_back(keys[rand() % keys.size()]);
    PlotPoints(locations);
    std::cout << "Calculating ..." << std::endl;
    auto results = TravellingTrojan(locations);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    CreateAnimation(results.second);
    if (results.second.size() != 0) {
      for (auto x : results.second[results.second.size()-1]) std::cout << x << std::endl;
      menu = "**************************************************************\n";
      std::cout << menu;
      std::cout << "The distance of the path is:" << results.first << std::endl;
      PlotPath(results.second[results.second.size()-1]);
    } else {
      std::cout << "The size of the path is 0" << std::endl;
    }
    menu = "**************************************************************\n"
           "You could find your animation at src/lib/output.avi.          \n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }
  case '5':
    break;
  default:
    std::cout << "Please select 1 - 5." << std::endl;
    PrintMenu();
    break;
  }
}


/**
 * CreateGraphFromCSVFile: Read the map data from the csv file
 * 
 */
void TrojanMap::CreateGraphFromCSVFile() {
  std::fstream fin;
  fin.open("src/lib/map.csv", std::ios::in);
  std::string line, word;

  getline(fin, line);
  while (getline(fin, line)) {
    std::stringstream s(line);

    Node n;
    int count = 0;
    while (getline(s, word, ',')) {
      word.erase(std::remove(word.begin(), word.end(), '\''), word.end());
      word.erase(std::remove(word.begin(), word.end(), '"'), word.end());
      word.erase(std::remove(word.begin(), word.end(), '['), word.end());
      word.erase(std::remove(word.begin(), word.end(), ']'), word.end());
      if (count == 0)
        n.id = word;
      else if (count == 1)
        n.lat = stod(word);
      else if (count == 2)
        n.lon = stod(word);
      else if (count == 3)
        n.name = word;
      else {
        word.erase(std::remove(word.begin(), word.end(), ' '), word.end());
        n.neighbors.push_back(word);
      }
      count++;
    }
    data[n.id] = n;
  }
  fin.close();
}

/**
 * PlotPoint: Given a location id, plot the point on the map
 * 
 * @param  {std::string} id : location id
 */
void TrojanMap::PlotPoint(std::string id) {
  std::string image_path = cv::samples::findFile("src/lib/input.jpg");
  cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
  auto result = GetPlotLocation(data[id].lat, data[id].lon);
  cv::circle(img, cv::Point(result.first, result.second), DOT_SIZE,
             cv::Scalar(0, 0, 255), cv::FILLED);
  cv::imshow("TrojanMap", img);
  cv::waitKey(1);
}
/**
 * PlotPoint: Given a lat and a lon, plot the point on the map
 * 
 * @param  {double} lat : latitude
 * @param  {double} lon : longitude
 */
void TrojanMap::PlotPoint(double lat, double lon) {
  std::string image_path = cv::samples::findFile("src/lib/input.jpg");
  cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
  auto result = GetPlotLocation(lat, lon);
  cv::circle(img, cv::Point(int(result.first), int(result.second)), DOT_SIZE,
             cv::Scalar(0, 0, 255), cv::FILLED);
  cv::startWindowThread();
  cv::imshow("TrojanMap", img);
  cv::waitKey(1);
}

/**
 * PlotPath: Given a vector of location ids draws the path (connects the points)
 * 
 * @param  {std::vector<std::string>} location_ids : path
 */
void TrojanMap::PlotPath(std::vector<std::string> &location_ids) {
  std::string image_path = cv::samples::findFile("src/lib/input.jpg");
  cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
  auto start = GetPlotLocation(data[location_ids[0]].lat, data[location_ids[0]].lon);
  cv::circle(img, cv::Point(int(start.first), int(start.second)), DOT_SIZE,
             cv::Scalar(0, 0, 255), cv::FILLED);
  for (auto i = 1; i < location_ids.size(); i++) {
    auto start = GetPlotLocation(data[location_ids[i - 1]].lat, data[location_ids[i - 1]].lon);
    auto end = GetPlotLocation(data[location_ids[i]].lat, data[location_ids[i]].lon);
    cv::circle(img, cv::Point(int(end.first), int(end.second)), DOT_SIZE,
               cv::Scalar(0, 0, 255), cv::FILLED);
    cv::line(img, cv::Point(int(start.first), int(start.second)),
             cv::Point(int(end.first), int(end.second)), cv::Scalar(0, 255, 0),
             LINE_WIDTH);
  }
  cv::startWindowThread();
  cv::imshow("TrojanMap", img);
  cv::waitKey(1);
}

/**
 * PlotPoints: Given a vector of location ids draws the points on the map (no path).
 * 
 * @param  {std::vector<std::string>} location_ids : points
 */
void TrojanMap::PlotPoints(std::vector<std::string> &location_ids) {
  std::string image_path = cv::samples::findFile("src/lib/input.jpg");
  cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
  for (auto x : location_ids) {
    auto result = GetPlotLocation(data[x].lat, data[x].lon);
    cv::circle(img, cv::Point(result.first, result.second), DOT_SIZE,
               cv::Scalar(0, 0, 255), cv::FILLED);
  }
  cv::imshow("TrojanMap", img);
  cv::waitKey(1);
}


/**
 * CreateAnimation: Create the videos of the progress to get the path
 * 
 * @param  {std::vector<std::vector<std::string>>} path_progress : the progress to get the path
 */
void TrojanMap::CreateAnimation(std::vector<std::vector<std::string>> path_progress){
  cv::VideoWriter video("src/lib/output.avi", cv::VideoWriter::fourcc('M','J','P','G'), 10, cv::Size(1248,992));
  for(auto location_ids: path_progress) {
    std::string image_path = cv::samples::findFile("src/lib/input.jpg");
    cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
    auto start = GetPlotLocation(data[location_ids[0]].lat, data[location_ids[0]].lon);
    cv::circle(img, cv::Point(int(start.first), int(start.second)), DOT_SIZE,
              cv::Scalar(0, 0, 255), cv::FILLED);
    for (auto i = 1; i < location_ids.size(); i++) {
      auto start = GetPlotLocation(data[location_ids[i - 1]].lat, data[location_ids[i - 1]].lon);
      auto end = GetPlotLocation(data[location_ids[i]].lat, data[location_ids[i]].lon);
      cv::circle(img, cv::Point(int(end.first), int(end.second)), DOT_SIZE,
                cv::Scalar(0, 0, 255), cv::FILLED);
      cv::line(img, cv::Point(int(start.first), int(start.second)),
              cv::Point(int(end.first), int(end.second)), cv::Scalar(0, 255, 0),
              LINE_WIDTH);
    }
    video.write(img);
    cv::imshow("TrojanMap", img);
    cv::waitKey(1);
  }
	video.release();
}
/**
 * GetPlotLocation: Transform the location to the position on the map
 * 
 * @param  {double} lat         : latitude 
 * @param  {double} lon         : longitude
 * @return {std::pair<double, double>}  : position on the map
 */
std::pair<double, double> TrojanMap::GetPlotLocation(double lat, double lon) {
  std::pair<double, double> bottomLeft(34.01001, -118.30000);
  std::pair<double, double> upperRight(34.03302, -118.26502);
  double h = upperRight.first - bottomLeft.first;
  double w = upperRight.second - bottomLeft.second;
  std::pair<double, double> result((lon - bottomLeft.second) / w * 1248,
                                   (1 - (lat - bottomLeft.first) / h) * 992);
  return result;
}

//-----------------------------------------------------
// TODO: Student should implement the following:
//-----------------------------------------------------
/**
 * GetLat: Get the latitude of a Node given its id.
 * 
 * @param  {std::string} id : location id
 * @return {double}         : latitude
 */
double TrojanMap::GetLat(std::string id) { 
  int lat = data[id].lat;
  return lat; 
}

/**
 * GetLon: Get the longitude of a Node given its id. 
 * 
 * @param  {std::string} id : location id
 * @return {double}         : longitude
 */
double TrojanMap::GetLon(std::string id) { 
  return data[id].lon;
}

/**
 * GetName: Get the name of a Node given its id.
 * 
 * @param  {std::string} id : location id
 * @return {std::string}    : name
 */
std::string TrojanMap::GetName(std::string id) { 
  return data[id].name;
}

/**
 * GetNeighborIDs: Get the neighbor ids of a Node.
 * 
 * @param  {std::string} id            : location id
 * @return {std::vector<std::string>}  : neighbor ids
 */
std::vector<std::string> TrojanMap::GetNeighborIDs(std::string id) {
    std::vector<std::string> result;
    result = data[id].neighbors;
    return result;
}


/**
 * CalculateDistance: Get the distance between 2 nodes. 
 * 
 * @param  {Node} a  : node a
 * @param  {Node} b  : node b
 * @return {double}  : distance in mile
 */
double TrojanMap::CalculateDistance(const Node &a, const Node &b) {
  // TODO: Use Haversine Formula:
  // dlon = lon2 - lon1;
  // dlat = lat2 - lat1;
  // a = (sin(dlat / 2)) ^ 2 + cos(lat1) * cos(lat2) * (sin(dlon / 2)) ^ 2;
  // c = 2 * arcsin(min(1, sqrt(a)));
  // distances = 3961 * c;
  double dlon = a.lon - b.lon;
  double dlat = a.lat - b.lat;
  double myaa = pow(sin(dlat/2),2)+cos(a.lat)*cos(b.lat)*pow(sin(dlon/2),2);
  double distance;
  if(1 > sqrt(myaa))
    distance = 3961*2*asin(sqrt(myaa));
  else
    distance = 3961*2*asin(sqrt(1));
  
  // where 3961 is the approximate radius of the earth at the latitude of
  // Washington, D.C., in miles
  return distance;
}

/**
 * CalculatePathLength: Calculates the total path length for the locations inside the vector.
 * 
 * @param  {std::vector<std::string>} path : path
 * @return {double}                        : path length
 */
double TrojanMap::CalculatePathLength(const std::vector<std::string> &path) {
  double sum = 0;
  int vsize = path.size();
  for(int i=1;i<vsize;i++){
    Node a = data[path[i-1]];
    Node b = data[path[i]];;
    double distance = CalculateDistance(a,b);
    sum += distance;
  }
  return sum;
}

/**
 * Autocomplete: Given a parital name return all the possible locations with
 * partial name as the prefix.
 *
 * @param  {std::string} name          : partial name
 * @return {std::vector<std::string>}  : a vector of full names
 */
std::vector<std::string> TrojanMap::Autocomplete(std::string name) {
  int nsize = name.size();
  std::map<std::string, Node>::iterator iter;
  std::vector<std::string> results;
  for(iter=data.begin();iter!=data.end();iter++){
    std::string nodename = iter->second.name;
    if(nodename.size() < nsize)
      continue;
    std::string tmp = nodename.substr(0,nsize);
    std::string tmpname = name;
    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
    std::transform(tmpname.begin(), tmpname.end(), tmpname.begin(), ::toupper);
    if(tmpname.compare(tmp)==0){
      results.insert(results.end(),nodename);
    }
  }
  return results;
}

std::string TrojanMap::GetID(std::string name){
  std::map<std::string, Node>::iterator iter;
  for(iter=data.begin();iter!=data.end();iter++){
    std::string nodename = iter->second.name;
    if(nodename.compare(name) == 0)
      return iter->second.id;
  }
  return "Not Found";
}

/**
 * GetPosition: Given a location name, return the position.
 *
 * @param  {std::string} name          : location name
 * @return {std::pair<double,double>}  : (lat, lon)
 */
std::pair<double, double> TrojanMap::GetPosition(std::string name) {
  std::map<std::string, Node>::iterator iter;
  std::pair<double, double> results(-1, -1);
  for(iter=data.begin();iter!=data.end();iter++){
    std::string curname = iter->second.name;
    if(curname.compare(name) == 0){
      results.first = iter->second.lat;
      results.second = iter->second.lon;
    }
  }
  return results;
}

// void dijk()

/**
 * CalculateShortestPath: Given 2 locations, return the shortest path which is a
 * list of id.
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 */
std::vector<std::string> TrojanMap::CalculateShortestPath(
    std::string location1_name, std::string location2_name) {
  std::string source_id = GetID(location1_name);
  std::string destin_id = GetID(location2_name);
  std::stack<std::string> s;
  int nsize = data.size();
  // Create a index map, from id to vector index
  std::map<std::string,int> indexlist;
  std::map<std::string,Node>::iterator iter;
  int i;
  for(i=0,iter=data.begin();iter!=data.end();iter++,i++){
    std::pair<std::string,int> newpair;
    newpair.first = iter->first;
    newpair.second = i;
    indexlist.insert(indexlist.end(),newpair);
  }
  // Create a visited distance vector and path frame
  std::vector<int> visited(nsize,0);
  std::vector<double> distance(nsize,INT_MAX);
  std::vector<std::vector<std::string>> path(nsize);
  // initial  
  distance[indexlist[source_id]] = 0;
  visited[indexlist[source_id]] = 1;
  path[indexlist[source_id]] = {source_id};
  int desindex = indexlist[destin_id];
  int srcindex = indexlist[source_id];
  s.push(source_id);
  // dijkstra algrithm
  while(visited[desindex] != 1){
    std::string tmpsid = s.top();
    s.pop();
    Node tmpnode = data[tmpsid];
    int tmpnsize = tmpnode.neighbors.size();
    for(int i=0;i<tmpnsize;i++){
      std::string tmpnid = tmpnode.neighbors[i];
      if(visited[indexlist[tmpnid]] == 1)
        continue;
      double tmpdistance = CalculateDistance(data[tmpsid],data[tmpnid]);
      if(distance[indexlist[tmpsid]]+tmpdistance < distance[indexlist[tmpnid]]){
        distance[indexlist[tmpnid]] = distance[indexlist[tmpsid]]+tmpdistance;
        path[indexlist[tmpnid]].clear();
        path[indexlist[tmpnid]] = path[indexlist[tmpsid]];
        path[indexlist[tmpnid]].insert(path[indexlist[tmpnid]].end(),tmpnid);
      }
    }
    double mindis = INT_MAX;
    int minindex;
    for(int i=0;i<nsize;i++){
      if(visited[i] == 1)
        continue;
      if(distance[i] < mindis ){
        mindis = distance[i];
        minindex = i;
      }
    }
    visited[minindex] = 1;
    for(std::map<std::string,int>::iterator it2 = indexlist.begin();it2!=indexlist.end();it2++) 
	  { 
		  if(it2->second==minindex){
			  s.push(it2->first);
        break;
      }
	  } 
  }
  std::vector<std::string> x;
  x = path[desindex];
  return x;
}

int fac(int i){
  if(i==0)
    return 1;
  else if(i<0)
    return 1;
  else
    return i*fac(i-1);
}

int Permutation_DFS(std::vector<int> mark, int step, int all, std::vector<std::vector<int>> &path,std::vector<int> &curpath)
{
  if (step == all){
    path.insert(path.end(),curpath);
    return 1;
  }
  int i;
  int total = 0;
  for (i = 1; i < all; i++)
  {
    if (mark[i])
      continue;
    mark[i] = 1;
    curpath[step] = i;
    total = total + Permutation_DFS(mark, step + 1, all, path,curpath);
    mark[i] = 0;
  }
  return total;
}

int isExist(std::vector<std::vector<int>> curpath,std::vector<int> allpath){
  int cpsize = curpath.size();
  return 0;
}

bool finishflag(std::vector<int> visited,std::vector<int> desindex){
  int dsize = desindex.size();
  for(int i=0;i<dsize;i++){
    if(visited[desindex[i]] == 1)
      continue;
    else
      return 0;
  }
  return 1;
}

std::vector<std::vector<std::string>> TrojanMap::CalculateShortestDistance(
    std::string location1_name, std::vector<std::string> location2_name,std::vector<double> &shortdis){
  std::string source_id = GetID(location1_name);
  int l2size = location2_name.size();
  std::vector<std::string> destin_id(l2size);
  for(int i=0;i<l2size;i++){
    destin_id[i] = GetID(location2_name[i]);
  }
  std::stack<std::string> s;
  int nsize = data.size();
  // Create a index map, from id to vector index
  std::map<std::string,int> indexlist;
  std::map<std::string,Node>::iterator iter;
  int i;
  for(i=0,iter=data.begin();iter!=data.end();iter++,i++){
    std::pair<std::string,int> newpair;
    newpair.first = iter->first;
    newpair.second = i;
    indexlist.insert(indexlist.end(),newpair);
  }
  // Create a visited distance vector and path frame
  std::vector<int> visited(nsize,0);
  std::vector<double> distance(nsize,INT_MAX);
  std::vector<std::vector<std::string>> path(nsize);
  // initial  
  distance[indexlist[source_id]] = 0;
  visited[indexlist[source_id]] = 1;
  path[indexlist[source_id]] = {source_id};
  std::vector<int> desindex(l2size);
  for(int i=0;i<l2size;i++){
    desindex[i] = indexlist[destin_id[i]];
  }
  // int desindex = indexlist[destin_id];
  int srcindex = indexlist[source_id];
  s.push(source_id);
  // dijkstra algrithm
  while(!finishflag(visited,desindex)){
    std::string tmpsid = s.top();
    s.pop();
    Node tmpnode = data[tmpsid];
    int tmpnsize = tmpnode.neighbors.size();
    for(int i=0;i<tmpnsize;i++){
      std::string tmpnid = tmpnode.neighbors[i];
      if(visited[indexlist[tmpnid]] == 1)
        continue;
      double tmpdistance = CalculateDistance(data[tmpsid],data[tmpnid]);
      if(distance[indexlist[tmpsid]]+tmpdistance < distance[indexlist[tmpnid]]){
        distance[indexlist[tmpnid]] = distance[indexlist[tmpsid]]+tmpdistance;
        path[indexlist[tmpnid]].clear();
        path[indexlist[tmpnid]] = path[indexlist[tmpsid]];
        path[indexlist[tmpnid]].insert(path[indexlist[tmpnid]].end(),tmpnid);
      }
    }
    double mindis = INT_MAX;
    int minindex;
    for(int i=0;i<nsize;i++){
      if(visited[i] == 1)
        continue;
      if(distance[i] < mindis ){
        mindis = distance[i];
        minindex = i;
      }
    }
    visited[minindex] = 1;
    for(std::map<std::string,int>::iterator it2 = indexlist.begin();it2!=indexlist.end();it2++) 
	  { 
		  if(it2->second==minindex){
			  s.push(it2->first);
        break;
      }
	  } 
  }
  std::vector<std::vector<std::string>> result;
  for(int i=0;i<l2size;i++){
    shortdis[i] = distance[desindex[i]];
    result.insert(result.end(),path[desindex[i]]);;
  }
  return result;
}
/**
 * Travelling salesman problem: Given a list of locations, return the shortest
 * path which visit all the places and back to the start point.
 *
 * @param  {std::vector<std::string>} input : a list of locations needs to visit
 * @return {std::pair<double, std::vector<std::vector<std::string>>} : a pair of total distance and the all the progress to get final path
 */
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan(
                                    std::vector<std::string> &location_ids) {
  std::pair<double, std::vector<std::vector<std::string>>> results;
  int lsize = location_ids.size();
  std::vector<std::vector<int>> path;
  int nsolu = fac(lsize-1)/2;
  // Create all solutions
  if(path.size() != nsolu){
    std::vector<int> mark(lsize,0);
    std::vector<std::vector<int>> tmppath;
    std::vector<int> curpath(lsize+1,0);
    Permutation_DFS(mark,1,lsize,tmppath,curpath);
    int psize = tmppath.size();
    // delete reverse solution
    for(int i=0;i<psize;i++){
      if(isExist(path,tmppath[i]) == 1)
        continue;
      path.insert(path.end(),tmppath[i]);
    }
  }
  // Create distance map
  std::vector<std::vector<double>> shortesdis;
  std::vector<std::vector<std::vector<std::string>>> pathlocid;
  // for(int i=0;i<lsize;i++){
  //   pathlocid[i].resize(lsize);
  // }
  for(int i=0;i<lsize;i++){
    std::vector<double> shortdis(lsize);
    std::string source_id = location_ids[i];
    std::vector<std::string> destin_id = location_ids;
    std::vector<std::vector<std::string>> results = CalculateShortestDistance(source_id,destin_id,shortdis);
    pathlocid.insert(pathlocid.end(),results);
    shortesdis.insert(shortesdis.end(),shortdis);
  }
  // Create 3D vector
  // std::vector<std::vector<std::vector<std::string>>> pathlocid;
  int psize = path.size();
  std::vector<double> sumdis(psize,0);
  std::vector<std::vector<std::string>> process(psize);
  for(int i=0;i<psize;i++){
    for(int j=1;j<lsize+1;j++){
      sumdis[i] += shortesdis[j-1][j];
    }
  }
  for(int i=0;i<psize;i++){
    for(int j=0;j<lsize+1;j++){
      process[i].insert(process[i].end(),location_ids[path[i][j]]);
    }
  }
  int minindex=0;
  for(int i=1;i<psize;i++){
    if(sumdis[minindex] > sumdis[i])
      minindex = i;
  }
  results.first = sumdis[minindex];
  results.second = process;
  location_ids = process[minindex];
  std::cout << results.second.size() << std::endl;
  return results;
} 

void TrojanMap::printindex() {
  std::map<std::string,int> indexlist;
  std::map<std::string,Node>::iterator iter;
  int i;
  for(i=0,iter=data.begin();iter!=data.end();iter++,i++){
    std::pair<std::string,int> newpair;
    newpair.first = iter->first;
    newpair.second = i;
    indexlist.insert(indexlist.end(),newpair);
  }
  std::map<std::string,int>::iterator it = indexlist.begin();
  // for(int time=0;time<100;time++){
  //   std::cout << it->first << " : " << it->second << std::endl;
  //   it++;
  // }
  for(it=indexlist.begin();it!=indexlist.end();it++){
    std::cout << it->first << " : " << it->second << std::endl; 
    iter++;
  }
}