/*
  geo1004.2023
  hw02 help code
  Hugo Ledoux <h.ledoux@tudelft.nl>
  2023-03-01
*/
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <typeinfo>
#include "definitions.h"
#include "geomtools.h"
#include <CGAL/Surface_mesh.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Simple_cartesian.h>


//-- https://github.com/nlohmann/json
//-- used to read and write (City)JSON
#include "json.hpp" //-- it is in the /include/ folder
using json = nlohmann::json;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_3<K> Tetrahedron;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::Surface_mesh<Point3> Surface_mesh;
typedef CGAL::Simple_cartesian<double>                           cartesian;
std::vector<Point3> get_coordinates(const json& j, bool translate = true);
void                save2obj(std::string filename, const json& j);
void                enrich_and_save(std::string filename, json& j);

using namespace std;

int main(int argc, const char * argv[]) {
    //-- will read the file passed as argument or 2b.city.json if nothing is passed
    const char* filename = (argc > 1) ? argv[1] : "../data/myfile.city.json";
    std::cout << "Processing: " << filename << std::endl;
    std::ifstream input(filename);
    json j;
    input >> j; //-- store the content of the file in a nlohmann::json object
    input.close();

    //-- convert each City Object in the file to OBJ and save to a file
    save2obj("out.obj", j);

    //-- enrich with some attributes and save to a new CityJSON
    enrich_and_save("out.city.json", j);

    return 0;
}


//-- write the OBJ file
void save2obj(std::string filename, const json& j) {
    std::ofstream ofile(filename);
    //-- fetch all the vertices in real-world coordinates (so "transform" is applied)
    std::vector<Point3> lspts = get_coordinates(j, true);
    for (auto& p : lspts) {
        ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }
    //-- iterate over each object in the file and output the CDT
    for (auto& co : j["CityObjects"].items()) {
        for (auto& g : co.value()["geometry"]) {
            if ( (g["type"] == "Solid") && (g["lod"] == "2.2") ) {   //-- LoD2.2 only!!!!!
                ofile << "o " << co.key() << std::endl;
                for (int i = 0; i < g["boundaries"].size(); i++) {
                    for (int j = 0; j < g["boundaries"][i].size(); j++) {
                        std::vector<std::vector<int>> gb = g["boundaries"][i][j];
                        std::vector<std::vector<int>> trs = construct_ct_one_face(gb, lspts);
                        for (auto& tr : trs) {
                            ofile << "f " << (tr[0] + 1) << " " << (tr[1] + 1) << " " << (tr[2] + 1) << std::endl;
                        }
                    }
                }
            }
        }
    }
    ofile.close();
    std::cout << "OBJ file written to disk: " << filename << std::endl;
}


/*
how we calculate volume:
step1 calcualte centrod, or anypoint in/on the boundary of polyhedradron, in case of any face/point is invalid, we use centroid for less effect from defects polyhetradron
notice* we need to assume that all the boudaries are triangulated and the orientation is correctly, aka. the normal vector point to the exterior(how to achieve it?
if the data input follow this rule)


deleted code:
CDT.insert(points.begin(),points.end());
          for (auto it = CDT.finite_cells_begin(); it != CDT.finite_cells_end(); ++it) {
          // Get the vertex indices of the current cell
            std::array<Point3, 4> vertex_indices;

            for (int i = 0; i < 4; ++i) {
            // Get the vertex handle of the i-th vertex of the cell
              auto vertex_handle = it->vertex(i);
              // Get the index of the vertex
              std::cout << "Vertex indices "<< i <<" f tetrahedron: " << vertex_handle->point() << "\n";
              vertex_indices[i] = vertex_handle->point();

            }
              // Do something with the vertex indices, e.g. print them out
              }



*/

// calculate single volume of tegrahedron
double determinant_4 (const double matrix[4][4], int n) {
    double det = 0;
    double submatrix[4][4];
    if (n == 2) return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
    else {
        for (int x = 0; x < n; x++) {
            int subi = 0;
            for (int i = 1; i < n; i++) {
                int subj = 0;
                for (int j = 0; j < n; j++) {
                    if (j == x) continue;
                    submatrix[subi][subj] = matrix[i][j];
                    subj++;
                }
                subi++;
            }
            int i = x;
            det = det + ((i % 2 == 0 ? 1 : -1) * matrix[0][x] * determinant_4( submatrix, n - 1 ));
        }
    }
    return det;
}
// sampling from surface
std::pair<Point3,vector<Point3>>  sampled_centroid( const double& total_area,const vector<double> area_weights,const std::vector<std::vector<int>> & all_triangulators,const std::vector<Point3>& lspts){

    typedef cartesian::Point_3             Point;
    int total_points = 2000;
    std::vector<Point3> sampled_points;
    std::vector<Point> points;
    
    for(int i = 0;i<all_triangulators.size();i++){
    std::vector<cartesian::Triangle_3> triangles;
    std::vector<int> tri = all_triangulators[i];
    // Create input triangles
    Point p1=Point(lspts[tri[0]].x(),lspts[tri[0]].y(),lspts[tri[0]].z());
    Point p2=Point(lspts[tri[1]].x(),lspts[tri[1]].y(),lspts[tri[1]].z());
    Point p3=Point(lspts[tri[2]].x(),lspts[tri[2]].y(),lspts[tri[2]].z());

    

    triangles.push_back(cartesian::Triangle_3(p1, p2, p3));
    CGAL::Random_points_in_triangles_3<Point> sample_g(triangles);
    int sample_nums = 2+int(total_points*area_weights[i]/total_area);
    std::copy_n(sample_g, sample_nums, std::back_inserter(points));
  
  // Create the generator, input is the vector of Triangle_3
  
    };
    
    // Get 100 random points in cdt
    cout<<"sample size"<<points.size()<<endl;
    Point3 centroid = Point3(0,0,0);
    for(Point pt:points){
        double x = pt.x();
        double y = pt.y();
        double z = pt.z();
        sampled_points.push_back(Point3(x,y,z));
        double centroidx = centroid.x();
        double centroidy = centroid.y();
        double centroidz = centroid.z();
        centroidx+=x/points.size();
        centroidy+=y/points.size();
        centroidz+=z/points.size();
        centroid = Point3(centroidx,centroidy,centroidz);


    }
    return std::make_pair(centroid,sampled_points);


}
double average_distance(const vector<Point3>& sampled_points,const Point3& centroid){
    double ptx = centroid.x();
    double pty = centroid.y();
    double ptz = centroid.z();
    double averate_dist = 0;
    for(auto& p1:sampled_points){

        double dx1;
        double dy1;
        double dz1;

        dx1 = p1.x()-ptx;
        dy1 = p1.y()-pty;
        dz1 = p1.z()-ptz;
        

        averate_dist+=std::sqrt(dx1*dx1+dy1*dy1+dz1*dz1);

    }
    return averate_dist/(sampled_points.size());
}

// not only volume, but also area
//calculate_tetrahedron(all_triangulators,vertices,centroid);
std::tuple<double,double,double> calculate_tetrahedron(const std::vector<std::vector<int>>& re,const std::vector<Point3>& lspts,const Point3& centroid){
    double centroidx = centroid.x();
    double centroidy = centroid.y();
    double centroidz = centroid.z();
    double total_volume=0;
    double total_area = 0;
    std::vector<double> area_weights;
    for(auto& tri: re){
        double matrix[4][4];
        Point3 p1=lspts[tri[0]];
        Point3 p2=lspts[tri[1]];
        Point3 p3=lspts[tri[2]];

        for(int pt=0;pt<3;pt++){
            Point3 vertex = lspts[tri[pt]];
            if(pt==0) p1 = vertex;
            if(pt==1) p2 = vertex;
            if(pt==2) p3 = vertex;
            //std::cout<<"vertex is "<<pt<<" "<<vertex<<std::endl;
            double ptx = vertex.x();
            double pty = vertex.y();
            double ptz = vertex.z();
            matrix[pt][0] =ptx;
            matrix[pt][1] =pty;
            matrix[pt][2] =ptz;
            matrix[pt][3] =1;
        }
        matrix[3][3] =1;
        matrix[3][2] =centroidz;
        matrix[3][1] =centroidy;
        matrix[3][0] =centroidx;



        /*for (auto& row : matrix) {
            for (auto& elem : row) {
                std::cout << elem << " ";
            }
            std::cout << std::endl;
        }*/
        double one_tri = determinant_4(matrix,4)/6;
        total_volume+=one_tri;




        double a = std::sqrt((p2.x() - p1.x()) * (p2.x() - p1.x()) + (p2.y() - p1.y()) * (p2.y() - p1.y()) + (p2.z() - p1.z()) * (p2.z() - p1.z()));
        double b = std::sqrt((p3.x() - p2.x()) * (p3.x() - p2.x()) + (p3.y() - p2.y()) * (p3.y() - p2.y()) + (p3.z() - p2.z()) * (p3.z() - p2.z()));
        double c = std::sqrt((p1.x() - p3.x()) * (p1.x() - p3.x()) + (p1.y() - p3.y()) * (p1.y() - p3.y()) + (p1.z() - p3.z()) * (p1.z() - p3.z()));

        // Calculate the semi-perimeter of the triangle
        double s = (a + b + c) / 2;

        // Calculate the area of the triangle using Heron's formula
        double area = std::sqrt(s * (s - a) * (s - b) * (s - c));
        total_area+=area;
        area_weights.push_back(area);


        //std::cout<<"volume of a tetrahedron is "<<one_tri<<std::endl;



    }
    Point3 centroid_new;
    std::vector<Point3> sampled_points;
    auto result = sampled_centroid(total_area,area_weights,re,lspts);
    centroid_new = result.first;
    sampled_points = result.second;
    std::cout<<"centroid is all tegrahedron is"<<centroid_new<<std::endl;
    double average_dist = average_distance(sampled_points,centroid_new);
    double roughIndex = 48.735*(average_dist*average_dist*average_dist)/(0.000001+total_volume+std::sqrt(total_area*total_area*total_area));

    return std::make_tuple(total_volume,total_area,roughIndex);
}
int calculate_aspect(const vector<int>& tri,const int& size_surface_semantics_array,const std::vector<Point3>& vertices){
    Point3 p1 = vertices[tri[0]];
    Point3 p2 = vertices[tri[1]];
    Point3 p3 = vertices[tri[2]];
    std::vector<double> vector1 = {p2.x()-p1.x(),p2.y()-p1.y(),p2.z()-p1.z()};
    std::vector<double> vector2 = {p3.x()-p1.x(),p3.y()-p1.y(),p3.z()-p1.z()};
    double norm_x = vector2[0]*vector1[0]+vector2[0]*vector1[1]+vector2[0]*vector1[2];
    double norm_y = vector2[1]*vector1[0]+vector2[1]*vector1[1]+vector2[1]*vector1[2];
    double length = sqrt(norm_x*norm_x+norm_y*norm_y);
    norm_y=norm_y/length;
    norm_x = norm_x/length;
    double slope = norm_y/norm_x;
    if(abs(norm_x-0)<=0.1 && abs(norm_y-0)<=0.1){
        return size_surface_semantics_array+8;

    }
    else if(norm_x==0 && norm_y>0){
        return size_surface_semantics_array+1;

    }
    else if(norm_x==0 && norm_y<0){
        return size_surface_semantics_array+5;

    }
    else if(slope>=1 && norm_y>0){
        return size_surface_semantics_array+1;

    }
    else if(slope>=1 && norm_y<0){
        return size_surface_semantics_array+5;

    }
    else if(slope<1 && slope>=0 && norm_x>0){
        return size_surface_semantics_array;

    }
    else if(slope<1 && slope>=0 && norm_x<0){
        return size_surface_semantics_array+4;

    }
    else if(slope<=-1 && norm_y>0){
        return size_surface_semantics_array+2;

    }
    else if(slope<=-1 && norm_y<0){
        return size_surface_semantics_array+6;

    }
    else if(slope<0 && slope>-1 && norm_x>0){
        return size_surface_semantics_array+7;

    }
    else if(slope<0 && slope>-1 && norm_x<0){
        return size_surface_semantics_array+3;

    }
    return size_surface_semantics_array+8;
}


std::tuple<double,double,double> calculate_volume_area(json& j,const json& co,const std::vector<Point3>& vertices){
    double volume_building = 0;
    double area_building = 0;
    double roughIndex = 0;
   
    double volume_building2=0,area_building2=0;
    if (co["type"] == "Building" && co["children"].size()>0) {
        // Print the type name
        //cout<<"run1"<<endl;
        std::vector<std::vector<int>> building_triangulators;
        //std::cout<<co["children"].size()<<std::endl;
        for(auto& child:co["children"]){

            //std::cout<<"geometry "<<j["CityObjects"][child]["geometry"][0]["boundaries"]<<std::endl;
            //std::cout<<"geometry type "<<j["CityObjects"][child]["geometry"][0]["type"]<<std::endl;
            //std::cout<<"geometry type "<<j["CityObjects"][child]["type"]<<std::endl;
            //std::cout<<"geometry type "<<j["CityObjects"][child]["semantics"]<<std::endl;
            //std::cout<<"geometry type "<<j["CityObjects"][child]["geometry"][0]["semantics"]<<std::endl;
            //std::cout<<"geometry boundaries size "<<j["CityObjects"][child]["geometry"][0]["boundaries"].size()<<std::endl;

            Point3 centroid(0,0,0);
            
            int number_of_point = 0;
            std::vector<std::vector<int>> all_triangulators;
            

            for(auto& solid:j["CityObjects"][child]["geometry"][0]["boundaries"]){
                
                
                json surface_semantics = j["CityObjects"][child]["geometry"][0]["semantics"];
                auto surfacesArray = surface_semantics["surfaces"];

                int roof_values = 0;
                cout<<"surface attributes are"<<j["CityObjects"][child]["geometry"][0]["semantics"]<<endl;
                for (auto it = surfacesArray.begin(); it != surfacesArray.end(); ) {
                            if (it.value()["type"] == "RoofSurface") {
                                roof_values = std::distance(surfacesArray.begin(), it);
                                it = surfacesArray.erase(it); 
                                break;// remove the object
                            } else {
                                     ++it;
                                    } }
                cout<<"roof_values "<<roof_values<<endl;
                int size_surface_semantics_array = surfacesArray.size();
                surfacesArray.push_back({ {"type", "RoofSurface"},{"orientation","EN" }}); //0
                surfacesArray.push_back({ { "type", "RoofSurface"},{"orientation", "NE" }}); //1
                surfacesArray.push_back({ { "type", "RoofSurface"},{"orientation", "NW" }}); //2
                surfacesArray.push_back({ { "type", "RoofSurface"},{"orientation", "WN" }}); //3
                surfacesArray.push_back({ { "type", "RoofSurface"},{"orientation", "WS" }}); //4
                surfacesArray.push_back({ { "type", "RoofSurface"},{"orientation", "SW" }}); //5
                surfacesArray.push_back({ { "type", "RoofSurface"},{"orientation", "SE" }}); //6
                surfacesArray.push_back({ { "type", "RoofSurface"},{"orientation", "ES" }}); //7
                surfacesArray.push_back({ { "type", "RoofSurface"},{"orientation", "Horizontal" }}); //8
                j["CityObjects"][child]["geometry"][0]["semantics"]["surfaces"]=surfacesArray;
                //cout<<j["CityObjects"][child]["geometry"][0]["semantics"]["surfaces"]<<endl;
                //cout<<j["CityObjects"][child]["geometry"][0]["semantics"]["values"]<<endl;
                auto& surfacesValues = j["CityObjects"][child]["geometry"][0]["semantics"]["values"][0];

                // start to iterate a cell on the boundaries
                int surface_index = 0;
                for(auto& surface:solid){
                    //std::cout<<"geometry boundaries size "<<surface<<std::endl;
                    std::vector<std::vector<int>> lsRings;
                    //cout<<"new surface is "<<surface<<endl;
                    //cout<<"surface size is "<<surface.size()<<endl;
                    for(auto& ring:surface){
                        cout<<surface<<endl;
                        std::vector<int> lsring;
                        //cout<<"new ring is "<<ring<<endl;
                        for(auto& pt:ring)
                        {   

                            //cout<<"new pt is "<<pt<<endl;
                            lsring.push_back(pt);
                            //cout<<"new pt been pushed is "<<pt<<endl;
                            Point3 pt3 = vertices[pt];
                            centroid=Point3(centroid.x() + pt3.x(), centroid.y() + pt3.y(),centroid.z() + pt3.z());
                            ++number_of_point;
                            //std::cout<<"point3 is "<<pt3<<std::endl;
                            //points.push_back(pt3);
                        }
                        //cout<<"new ring been pushed is "<<ring<<endl;
                        lsRings.push_back(lsring);
                        //cout<<"new ring is"<<endl;
                    }
                    //cout<<"current Rings is"<<endl;
                    std::vector<std::vector<int>> constrain_triangulors = construct_ct_one_face( lsRings,  vertices);
                    //cout<<"current Rings is"<<lsRings.size()<<endl;
                    // add to all_triangulators volume
                    building_triangulators.insert(building_triangulators.end(), constrain_triangulors.begin(), constrain_triangulors.end());
                    //cout<<"current triangulars is"<<constrain_triangulors.size()<<endl;
                    //if(constrain_triangulors.size() == 0){constrain_triangulors = lsRings;}
                    all_triangulators.insert(all_triangulators.end(), constrain_triangulors.begin(), constrain_triangulors.end());
                    int surface_value =surfacesValues[surface_index];
                    
                    if(surface_value==roof_values){
                        //cout<<"roof values"<<surface_value<<endl;
                        //cout<<"run5"<<endl;
                        int new_value= 0;
                        if(constrain_triangulors.size()==0){int new_value = calculate_aspect(lsRings[0],size_surface_semantics_array,vertices);}
                        else{
                            new_value = calculate_aspect(constrain_triangulors[0],size_surface_semantics_array,vertices);
                        }

                        
                        surfacesValues[surface_index] = new_value;
                        cout<<"new_value is"<<new_value<<endl;
                        //cout<<"next step"<<endl;
                        

                    }
                    //cout<<"run2"<<endl;
                    //cout<<j["CityObjects"][child]<<endl;
                    //cout<<"new value"<<surface_index<<endl;
                    surface_index+=1;
                    //cout<<"new value"<<surface<<endl;
                    //cout<<"solid is"<<solid<<endl;
                }
                
            }
            //cout<<"run3"<<endl;
            //get_centroid_point
             centroid=Point3(centroid.x()/number_of_point, centroid.y()/number_of_point,centroid.z()/number_of_point);
                //std::cout<<"centroid is "<<centroid<<std::endl;
                // the total volume
                // in fact we don't need centroid
            
            auto result = calculate_tetrahedron(all_triangulators,vertices,centroid);
                //auto result2 = calculate_tetrahedron(all_triangulators,vertices,pt_zero);
            
            volume_building += std::get<0>(result);
            area_building += std::get<1>(result);
            
        }
        Point3 pt_zero(0,0,0);
        auto result2 = calculate_tetrahedron(building_triangulators,vertices,pt_zero);
        //std::cout<<"building seperate equal? "<<volume_building<<"  "<<(std::get<0>(result2))<<std::endl;
        //std::cout<<"area centorid equal? "<<area_building<<"  "<<(std::get<1>(result2))<<std::endl;
        
        area_building2 = std::get<1>(result2);
        volume_building2 = std::get<0>(result2);

        roughIndex = std::get<2>(result2);
        //break;
        //co["attributes"]["volume"] = dis(gen);
    }
    return std::make_tuple(volume_building,area_building,roughIndex);



}

double calculate_OBB_Volumn(const json& j,const json& co,const std::vector<Point3>& vertices)
{
    double volume_building = 0;
    double area_building = 0;
    double b=0;
    vector<Point3> BBD_pt;
    if (co["type"] == "Building") {
        // Print the type name
        //std::cout<<co["children"].size()<<std::endl;
        for(auto& child:co["children"]){
            int number_of_point = 0;
            vector<int> lsRings;
            for(auto& solid:j["CityObjects"][child]["geometry"][0]["boundaries"]){
                // start to iterate a cell on the boundaries
                for(auto& surface:solid){
                    //std::cout<<"geometry boundaries size "<<surface<<std::endl;
                    for(auto& ring:surface){
                        for(auto& pt:ring)
                        {
                            lsRings.push_back(pt);
                            ++number_of_point;
                            //std::cout<<"point3 is "<<pt3<<std::endl;
                            //points.push_back(pt3);
                        }
                    }}}

            for (auto &tri:lsRings) {
                Point3 vertex_temp = vertices[tri];
                BBD_pt.push_back(vertex_temp);
            }
        }
        std::array<Point3, 8> obb_points;
        CGAL::oriented_bounding_box(BBD_pt, obb_points);
//        for (int i=0;i<8;i++){
//            cout<<"第"<<i<<" 个点坐标为"<<obb_points[i];
//        }
        Surface_mesh obb_sm;
        CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                              obb_points[4], obb_points[5], obb_points[6], obb_points[7], obb_sm);
        PMP::triangulate_faces(obb_sm);
        auto v_obb=PMP::volume(obb_sm);
        b += static_cast<double>(v_obb);
        std::cout << "Volumn of OBB: " << b << std::endl;
    }
return b;
}


Point3  calculate_centriod(const json& j,const json& co,const std::vector<Point3>& vertices){
    double volume_building = 0;
    double area_building = 0;Point3 centriod;

        // Print the type name
        vector<Point3>sample_points;
        //std::cout<<co["children"].size()<<std::endl;
        double avgarea;double maxarea=0;double minarea=99999;int num_tri=0;

        for(auto& child:co["children"]){

            //std::cout<<"geometry "<<j["CityObjects"][child]["geometry"][0]["boundaries"]<<std::endl;
            //std::cout<<"geometry type "<<j["CityObjects"][child]["geometry"][0]["type"]<<std::endl;
            //std::cout<<"geometry type "<<j["CityObjects"][child]["type"]<<std::endl;
            //std::cout<<"geometry type "<<j["CityObjects"][child]["semantics"]<<std::endl;
            //std::cout<<"geometry type "<<j["CityObjects"][child]["geometry"][0]["semantics"]<<std::endl;
            std::cout<<"geometry boundaries size "<<j["CityObjects"][child]["geometry"][0]["boundaries"].size()<<std::endl;
            Point3 centroid(0,0,0);
            std::vector<std::vector<int>> all_triangulators;
            int number_of_point = 0;
            for(auto& solid:j["CityObjects"][child]["geometry"][0]["boundaries"]){
                // start to iterate a cell on the boundaries
                for(auto& surface:solid){
                    //std::cout<<"geometry boundaries size "<<surface<<std::endl;
                    std::vector<std::vector<int>> lsRings;
                    for(auto& ring:surface){
                        std::vector<int> lsring;
                        for(auto& pt:ring)
                        {
                            lsring.push_back(pt);
                            Point3 pt3 = vertices[pt];
                            centroid=Point3(centroid.x() + pt3.x(), centroid.y() + pt3.y(),centroid.z() + pt3.z());
                            ++number_of_point;
                            //std::cout<<"point3 is "<<pt3<<std::endl;
                            //points.push_back(pt3);
                        }
                        lsRings.push_back(lsring);
                    }
                    std::vector<std::vector<int>> constrain_triangulors = construct_ct_one_face( lsRings,  vertices);
                    // add to all_triangulators volume
                    all_triangulators.insert(all_triangulators.end(), constrain_triangulors.begin(), constrain_triangulors.end());
                }
            }


            for(auto& tri: all_triangulators){
                Point3 p1=vertices[tri[0]];Point3 p2=vertices[tri[1]];Point3 p3=vertices[tri[2]];
                double a = std::sqrt((p2.x() - p1.x()) * (p2.x() - p1.x()) + (p2.y() - p1.y()) * (p2.y() - p1.y()) + (p2.z() - p1.z()) * (p2.z() - p1.z()));
                double b = std::sqrt((p3.x() - p2.x()) * (p3.x() - p2.x()) + (p3.y() - p2.y()) * (p3.y() - p2.y()) + (p3.z() - p2.z()) * (p3.z() - p2.z()));
                double c = std::sqrt((p1.x() - p3.x()) * (p1.x() - p3.x()) + (p1.y() - p3.y()) * (p1.y() - p3.y()) + (p1.z() - p3.z()) * (p1.z() - p3.z()));

                // Calculate the semi-perimeter of the triangle
                double s = (a + b + c) / 2;

                // Calculate the area of the triangle using Heron's formula
                double tri_area = std::sqrt(s * (s - a) * (s - b) * (s - c));
                avgarea+=tri_area;
                if (tri_area>maxarea)
                    maxarea=tri_area;
                if(tri_area<minarea)
                    minarea=tri_area;
                num_tri++;
            }
        }
        avgarea=avgarea/num_tri;

        for(auto& child:co["children"]){
            std::vector<std::vector<int>> all_triangulators;
            int number_of_point = 0;
            for(auto& solid:j["CityObjects"][child]["geometry"][0]["boundaries"]){
                // start to iterate a cell on the boundaries
                for(auto& surface:solid){
                    //std::cout<<"geometry boundaries size "<<surface<<std::endl;
                    std::vector<std::vector<int>> lsRings;
                    for(auto& ring:surface){
                        std::vector<int> lsring;
                        for(auto& pt:ring)
                        {
                            lsring.push_back(pt);
                            Point3 pt3 = vertices[pt];
                            ++number_of_point;
                            //std::cout<<"point3 is "<<pt3<<std::endl;
                            //points.push_back(pt3);
                        }
                        lsRings.push_back(lsring);
                    }
                    std::vector<std::vector<int>> constrain_triangulors = construct_ct_one_face( lsRings,  vertices);
                    // add to all_triangulators volume
                    all_triangulators.insert(all_triangulators.end(), constrain_triangulors.begin(), constrain_triangulors.end());
                }
            }
            for(auto& tri: all_triangulators){
                Point3 p1=vertices[tri[0]];Point3 p2=vertices[tri[1]];Point3 p3=vertices[tri[2]];
                double a = std::sqrt((p2.x() - p1.x()) * (p2.x() - p1.x()) + (p2.y() - p1.y()) * (p2.y() - p1.y()) + (p2.z() - p1.z()) * (p2.z() - p1.z()));
                double b = std::sqrt((p3.x() - p2.x()) * (p3.x() - p2.x()) + (p3.y() - p2.y()) * (p3.y() - p2.y()) + (p3.z() - p2.z()) * (p3.z() - p2.z()));
                double c = std::sqrt((p1.x() - p3.x()) * (p1.x() - p3.x()) + (p1.y() - p3.y()) * (p1.y() - p3.y()) + (p1.z() - p3.z()) * (p1.z() - p3.z()));
                //cout<<p1<<p2<<p3<<endl;
                // Calculate the semi-perimeter of the triangle
                double s = (a + b + c) / 2;

                // Calculate the area of the triangle using Heron's formula
                double tri_area = std::sqrt(s * (s - a) * (s - b) * (s - c));
                int p_num;
                if (tri_area>avgarea)
                    p_num=4;
                if(tri_area<avgarea)
                    p_num=2;
                //cout<<"the current triangle area is "<<tri_area<<endl;

                    // Create a random number generator
                CGAL::Random rng;
                // Create a point generator
                CGAL::Random_points_in_triangle_3<Point3> generator(p1, p2, p3, rng);
                // Generate 10 random points within the triangle
                for (int i = 0; i < p_num; i++) {
                    Point3 p = *generator++;
                    sample_points.push_back(p);
                    //std::cout << p << std::endl;
                }
            }
        }

        //break;
        //co["attributes"]["volume"] = dis(gen);
        int num_pt=0;
        double avgx=0;double avgy=0;double avgz=0;
        cout<<sample_points.size();
        for (auto&pt:sample_points){
            //cout<<pt<<endl;
            avgx+=pt.x();avgy+=pt.y();avgz+=pt.z();
            num_pt++;
        }
        cout<<"Point number is: "<<num_pt<<endl;
        centriod=Point3(avgx/num_pt,avgy/num_pt,avgz/num_pt);

    return centriod;
}


//-- add a new attribute "volume" to each City Object and assign a random value
void enrich_and_save(std::string filename, json& j) {
    //-- seed to generate a random number
    //-- https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(100.0, 200.0);
    //-- add an attribute "volume"
    std::vector<Point3> vertices = get_coordinates(j,  1);
    for (auto& co : j["CityObjects"]) {
        if (co["type"] == "Building"){
            auto result = calculate_volume_area(j, co,vertices);
            double volume_building = std::get<0>(result);
            double area_building = std::get<1>(result);
            double roughIndex = std::get<2>(result);
            auto temp=calculate_OBB_Volumn(j, co,vertices);
            double rectan=volume_building/temp;
            std::cout<<"the building volume is "<<volume_building<<std::endl;
            std::cout<<"the oreintBB volume is "<<temp<<endl;
            co["attributes"]["volume"] = volume_building;
            co["attributes"]["roughIndex"] =roughIndex;
            co["attributes"]["rectangularity"] = rectan;
            cout<<"the rectangularity is "<<rectan<<endl;
            co["attributes"]["area"] = area_building;
            double hemisphe=3*sqrt(2*M_PI)*volume_building/(sqrt(area_building*area_building*area_building));
            cout<<"the hemisphericality is "<<hemisphe<<endl;
            co["attributes"]["hemisphericality"] = hemisphe;
//            Point3 centriod=calculate_centriod(j, co,vertices);
//            cout<<"cebtriod Point is "<<centriod<<endl;

        }
    }
    //-- write to disk the modified city model (myfile.city.json)
    std::ofstream o(filename);
    o << j.dump(2) << std::endl;
    o.close();
    std::cout << "Enriched CityJSON file written to disk: " << filename << std::endl;
}


//-- get real-world coordinates from the "vertices" property
//-- https://www.cityjson.org/specs/#transform-object
//-- param translate is to use the translation in the "transform",
//-- it can be put to false to make the coords smaller (and better for computations)
std::vector<Point3> get_coordinates(const json& j, bool translate) {
    std::vector<Point3> lspts;    std::vector<std::vector<int>> lvertices = j["vertices"];
    if (translate) {
        for (auto& vi : lvertices) {
            double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
            double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
            double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
            lspts.push_back(Point3(x, y, z));
        }
    } else {
        //-- do not translate, useful to keep the values low for downstream processing of data
        for (auto& vi : lvertices) {
            double x = (vi[0] * j["transform"]["scale"][0].get<double>());
            double y = (vi[1] * j["transform"]["scale"][1].get<double>());
            double z = (vi[2] * j["transform"]["scale"][2].get<double>());
            lspts.push_back(Point3(x, y, z));
        }
    }
    return lspts;
}