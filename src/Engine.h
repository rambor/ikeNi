//
// Created by xos81802 on 07/09/2021.
//
// Copyright 2021 Robert P. Rambo
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// 
//

#ifndef MIRUSAXS_ENGINE_H
#define MIRUSAXS_ENGINE_H


#include <string>
#include <cmath>
#include <sastools/include/PDBModel.h>
#include <utility>
#include <cstdlib>

class Engine {

    struct Pose {
        float theta, phi, psi, chi;
        vector3 vec;
        Pose(float t, float p, float ps, float ch, vector3 vec) : theta(t), phi(p), psi(ps), chi(ch), vec(vec) {
        }
    };


    std::string reference_filename, target_filename;
    PDBModel reference_model;//, target_model;

    std::vector<vector3> centered_ref_coordinates;

    void extractCoordinates();

public:
    Engine(){std::cout << "empty -->" << std::endl;};
    Engine(std::string reference, std::string target_filename);
    ~Engine()= default;

    Engine(const Engine & temp){
        reference_filename = temp.reference_filename;
        target_filename = temp.target_filename;

        for(auto & vec : temp.centered_ref_coordinates){
            centered_ref_coordinates.emplace_back(vector3(vec));
        }

        reference_model = std::move(reference_model);
    }


    void search();

    float getEpsilon(int total);

    void writeCoordinatesToFile(std::vector<vector3> &vecs, std::string name);

    void generateRotationMatrix(float theta, float phi, float psi, std::vector<vector3> & rotationMatrix);

    void rotateTarget(std::vector<vector3> & rotationMatrix, std::vector<vector3> & target, std::vector<vector3> & posed);

    void translateTarget(vector3 &translation_vector, std::vector<vector3> &target);

    float getChiFromFile(std::string filename);

    float closestDistance(std::vector<vector3> &pose);

    void writeMergedCoordinates(std::string name, PDBModel &targetModel, std::vector<vector3> & posed);
};


#endif //MIRUSAXS_ENGINE_H
