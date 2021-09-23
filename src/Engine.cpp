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


#include "Engine.h"



Engine::Engine(std::string reference, std::string target_filename) : reference_filename(std::move(reference)), target_filename(std::move(target_filename)) {

}

void Engine::search(){

    reference_model = PDBModel(reference_filename, false, true);
    reference_model.setCVXHullPoints();
    reference_model.writeCVXPointsCoordinatesToFile();
    reference_model.writeCenteredCoordinatesToFile("centered_reference");

    PDBModel target_model = PDBModel(target_filename, false, true);
    unsigned int total_in_target = target_model.getTotalCoordinates();

    unsigned int totalCoords = reference_model.getTotalCoordinates();

    const float * cX = reference_model.getCenteredXVec();
    const float * cY = reference_model.getCenteredYVec();
    const float * cZ = reference_model.getCenteredZVec();

    for(unsigned int i = 0; i<totalCoords; i++){
        //std::cout << reference_model.getCenteredXVec()[i] << " - " << reference_model.getCenteredYVec()[i] << " " << reference_model.getCenteredZVec()[i] << std::endl;
        //    if (reference_model.belongsToResidue(i) && reference_model.isBackbone(i)){
        centered_ref_coordinates.emplace_back(
                vector3(cX[i], cY[i], cZ[i])
        );
        //    }
    }

    // create rotation points on surface of sphere
    std::cout << "        SEARCH :: SPHERICAL GRID"  << std::endl;
    //std::cout << "  total trials :: " << (theta_limit*phi_limit*psi_limit) << std::endl;

    totalCoords = centered_ref_coordinates.size();
    int total_sphere_grid_points = 131;
    //std::cout << "  total trials :: " << (total_sphere_grid_points*psi_limit) << std::endl;
    // grid search
    // use Fibonacci grid
    auto golden_ratio = (1.0f + sqrtf(5.0))/2.0;
    auto inv_golden_ratio = (float)((2.0f*M_PI)/golden_ratio);

    float epsilon = getEpsilon(total_sphere_grid_points), projection, perp, magnitude, inv_sq_length, dot;
    //auto delta_psi = (float)(2.0*M_PI/(double)psi_limit);
    int totalinGrid = 0;
    const vector3 * vec;
    vector3 a1, a2;

    // target
    float target_radius = target_model.getSMax();
    float adjusted_magnitude;

    unsigned int theta_limit = 37;
    unsigned int phi_limit = 19;
    unsigned int psi_limit = 17;
    auto delta_theta = (float)(2.0*M_PI/(double)theta_limit);
    auto delta_phi = (float)(2.0*M_PI/(double)phi_limit);
    auto delta_psi = (float)(2.0*M_PI/(double)psi_limit);

    Pose best_pose = Pose(0,0,0,0, vector3(0,0,0));

    std::vector<vector3> rotationMatrix(3);
    std::vector<vector3> rotation_centers;
    std::vector<vector3> fibonacci_lines;
    std::vector<int> kept_indices_in_reference;
    std::vector<vector3> posed_target(target_model.getTotalCoordinates());
    std::vector<vector3> copy_posed_target(target_model.getTotalCoordinates());
    target_model.setVectorModelToCenter();
    auto target_coord_vecs = *(target_model.getModelVector());
    target_model.writeTranslatedCoordinatesToFile("centered_target", *target_model.getModelVector());

    int model_count=0;
    float temp_chi, current_chi = FLT_MAX, dis;

    // precompute spherical grid vectors
    for(int j=0; j < total_sphere_grid_points; j++) {

        float angle_theta = j * inv_golden_ratio; // rotation around X
        float angle_phi = acosf(
                1.0f - 2.0f * (j + epsilon) / (total_sphere_grid_points - 1 + 2.0f * epsilon)); // rotation around z

        // want to rotate about the vector
        // find the point in reference model that is closest to line and furthest away from center
        fibonacci_lines.emplace_back(vector3(sinf(angle_phi) * cosf(angle_theta), sinf(angle_phi) * sinf(angle_theta),
                                             cosf(angle_phi)));

        vector3 line = fibonacci_lines[j];

        inv_sq_length = 1.0f / line.sqlength();

        float max_proj = 0.0f;
        float min_perp = 1000.0f;
        int kept_index = 0;
        // find coordinate that is closest to project Fibonacci line along grid point
        for (int c = 0; c < totalCoords; c++) {
            vec = &centered_ref_coordinates[c];
            dot = line.dot(*vec);
            a1 = line * (dot * inv_sq_length);
            projection = a1.length();
            a2 = *vec - a1;
            perp = a2.length();
            if (projection > max_proj && perp < min_perp) {
                max_proj = projection;
                min_perp = perp;
                kept_index = c;
            }
        }
        kept_indices_in_reference.push_back(kept_index);
    }

    std::string filename, crysol_file = "crysol_summary.txt";
    std::string kept_name = "summary_kept.txt";
    FILE * pFile;
    pFile = fopen(kept_name.c_str(), "w");
    fprintf(pFile, "REMARK model_count theta phi psi CHI2 trans_vector\n");
    fclose(pFile);

    // for each grid point - do random sampling of orientations


    for(int j=0; j < total_sphere_grid_points; j++){
        // for each grid point try 10 different poses?

    }
    // score model -> j index

    // first pass try each


    // update probability models for theta, psi, and phi (sigma and mu)
    // each j_index has its own theta, psi, phi model



//    for(int j=0; j < total_sphere_grid_points; j++){
//
//        float angle_theta = j*inv_golden_ratio; // rotation around X
//        float angle_phi = acosf(1.0f - 2.0f*(j+epsilon)/(total_sphere_grid_points - 1 + 2.0f*epsilon)); // rotation around z
//
//        // want to rotate about the vector
//        // find the point in reference model that is closest to line and furthest away from center
//        vector3 line = vector3(sinf(angle_phi)*cosf(angle_theta), sinf(angle_phi)*sinf(angle_theta), cosf(angle_phi));
//
//        inv_sq_length = 1.0f/line.sqlength();
//
//        float max_proj = 0.0f;
//        float min_perp = 1000.0f;
//        int kept_index = 0;
//        // find coordinate that is closest to project Fibonacci line along grid point
//        for (int c=0; c<totalCoords; c++){
//            vec = &centered_ref_coordinates[c];
//            dot = line.dot(*vec);
//            a1 = line*(dot*inv_sq_length);
//            projection = a1.length();
//            a2 = *vec - a1;
//            perp = a2.length();
//            if (projection > max_proj && perp < min_perp){
//                max_proj = projection;
//                min_perp = perp;
//                kept_index = c;
//            }
//        }
//
//        // get vector of ref atom closest to the projected fibonacci grid point
//        vec = &centered_ref_coordinates[kept_index];

        // take the target and perform rotations
        for (unsigned int theta_n=0; theta_n<theta_limit; theta_n++){ // range is from 0 to 180
            float angle_theta = theta_n*delta_theta;
            for (unsigned int phi_n=0; phi_n<phi_limit; phi_n++){ // range is from 0 to 360
                float angle_phi = phi_n*delta_phi;
                for (unsigned int psi_n=0; psi_n<psi_limit; psi_n++){ // range is from 0 to 360
                    generateRotationMatrix(angle_theta, angle_phi, psi_n*delta_psi, rotationMatrix);
                    //rotate target
                    rotateTarget(rotationMatrix, target_coord_vecs, copy_posed_target);

                    // translate to new position
                    for(int j=0; j < total_sphere_grid_points; j++){
                        vec = &centered_ref_coordinates[kept_indices_in_reference[j]];
                        adjusted_magnitude = vec->length() + target_radius;
                        vector3 line = fibonacci_lines[j];
                        a1 = line.normalize()*adjusted_magnitude;

                        std::copy(copy_posed_target.begin(), copy_posed_target.end(), posed_target.begin());

                        translateTarget(a1, posed_target);
                        // test if too far or too close
                        dis=closestDistance(posed_target);

                        if (dis < 1){
                            while (dis < 1){
                                std::cout << "Too close, adjusting " << dis << std::endl;
                                a1 = line.normalize()*2*dis;
                                translateTarget(a1, posed_target);
                                dis=closestDistance(posed_target);
                            }
                        } else if (dis > 3){
                            while (dis > 3){
                                a1 = line.normalize()*(-0.5*dis); // sub
                                std::cout << "Too far, adjusting " << a1.length() << " " << std::endl;
                                translateTarget(a1, posed_target);
                                dis=closestDistance(posed_target);
                            }
                        }

                        // rotate the target
                        //filename = "pose_" + std::to_string(model_count);
                        //target_model.writeTranslatedCoordinatesToFile("pose", posed_target);
                        writeMergedCoordinates("pose", target_model, posed_target);

                        // test pose
                        std::system("crysol pose.pdb Complex_Moore_sx.dat");
                        std::system("rm pose00.fit");
                        std::system("rm pose00.log");
                        // get chi value
                        temp_chi = getChiFromFile(crysol_file);

                        if (temp_chi < 2 || temp_chi < current_chi){

                            if(temp_chi < current_chi){
                                std::system("cp pose.pdb best_pose.pdb");
                            }

                            current_chi = temp_chi;
                            best_pose.theta = angle_theta;
                            best_pose.phi = angle_phi;
                            best_pose.psi = psi_n*delta_psi;
                            best_pose.chi = current_chi;
                            best_pose.vec = line;
                            std::string pose_name = "kept_"+std::to_string(model_count) + ".pdb";
                            std::string commaand = "mv pose.pdb " + pose_name;
                            std::system(commaand.c_str());

                            pFile = fopen(kept_name.c_str(), "a");
                            fprintf(pFile, "%8d %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                                    model_count,
                                    best_pose.theta,
                                    best_pose.phi,
                                    best_pose.psi,
                                    best_pose.vec.x,
                                    best_pose.vec.y,
                                    best_pose.vec.z,
                                    best_pose.chi
                            );
                            fclose(pFile);
                        }

                        model_count++;
                    }

//                    adjusted_magnitude = vec->length() + target_radius;
//                    a1 = line.normalize()*adjusted_magnitude;
//                    //rotation_centers.emplace_back(vector3(a1));
//                    translateTarget(a1, posed_target);
//                    // test if too far or too close
//                    dis=closestDistance(posed_target);
//
//                    if (dis < 1){
//                        while (dis < 1){
//                            std::cout << "Too close, adjusting " << dis << std::endl;
//                            a1 = line.normalize()*2*dis;
//                            translateTarget(a1, posed_target);
//                            dis=closestDistance(posed_target);
//                        }
//                    } else if (dis > 3){
//                        while (dis > 3){
//                            a1 = line.normalize()*(-0.5*dis); // sub
//                            std::cout << "Too far, adjusting " << a1.length() << " " << std::endl;
//                            translateTarget(a1, posed_target);
//                            dis=closestDistance(posed_target);
//                        }
//                    }
//
//                    // rotate the target
//                    //filename = "pose_" + std::to_string(model_count);
//                    //target_model.writeTranslatedCoordinatesToFile("pose", posed_target);
//                    writeMergedCoordinates("pose", target_model, posed_target);
//
//                    // test pose
//                    std::system("crysol pose.pdb Complex_Moore_sx.dat");
//                    std::system("rm pose00.fit");
//                    std::system("rm pose00.log");
//                    // get chi value
//                    temp_chi = getChiFromFile(crysol_file);
//                    if (temp_chi < current_chi){
//                        current_chi = temp_chi;
//                        best_pose.theta = angle_theta;
//                        best_pose.phi = angle_phi;
//                        best_pose.psi = psi_n*delta_psi;
//                        best_pose.chi = current_chi;
//                        best_pose.trans = a1.length();
//                        std::system("mv pose.pdb kept.pdb");
//                        pFile = fopen(kept_name.c_str(), "a");
//                        fprintf(pFile, "%8d %5.2f %5.2f %5.2f %5.2f %.1f\n", model_count, best_pose.theta, best_pose.phi, best_pose.psi, best_pose.chi, best_pose.trans);
//                        fclose(pFile);
//                    }
//                    model_count++;
                }
                //std::cout << "Theta " << angle_theta << " PHI " << angle_phi << std::endl;
            }
//        }

    }

//    writeCoordinatesToFile(rotation_centers, "centers");

//std::cout << "--" <<std::endl;
//    for (int c=0; c<totalCoords; c++){
//        vec = &centered_ref_coordinates[c];
//        std::cout << vec->x << " " << vec->y << " "  << vec->z << std::endl;
//    }


}


float Engine::getEpsilon(int total){

    if (total >= 600000){
        return 214;
    } else if (total >= 400000){
        return 75;
    } else if (total >= 11000){
        return 27;
    } else if (total >= 890){
        return 10;
    } else if (total >= 177){
        return 3.33;
    } else if (total >= 24){
        return 1.33;
    }

    return 0.33;
}


void Engine::writeCoordinatesToFile(std::vector<vector3> & vecs, std::string name) {

    std::string residue_index;
    name = name + ".pdb";
    //const char * outputFileName = name.c_str() ;
    //const char * originalPDBFilename = this->filename.c_str();
    FILE * pFile;
    pFile = fopen(name.c_str(), "w");
    unsigned int totalAtoms = vecs.size();
    vector3 * vec;

    fprintf(pFile,"REMARK  CENTERED COORDINATES : %s\n", "Perimeter");
    for (unsigned int n=0; n < totalAtoms; n++) {
        vec = &vecs[n];
        residue_index = boost::lexical_cast<std::string>(n+1);
        //fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", n+1, trimmedAtomType[n].c_str(), resi[n].c_str(), chainID[n].c_str(), residue_index.c_str(), centeredX[n], centeredY[n], centeredZ[n] );
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", n+1, "CA", "ALA", "A", residue_index.c_str(), vec->x, vec->y, vec->z );
    }
    fprintf(pFile,"END\n");
    fclose(pFile);
}


void Engine::generateRotationMatrix(float theta, float phi, float psi, std::vector<vector3> & rotationMatrix){

    //float twoPI = (float)(2.0*M_PI);
    float costheta = cos(theta);
    float sintheta = sin(theta);
    float cosphi = cos(phi);
    float sinphi = sin(phi);
    float cospsi = cos(psi);
    float sinpsi = sin(psi);

    float sinphi_costheta= sinphi*costheta;
    float sinphi_sintheta= sinphi*sintheta;

//    Eigen::Matrix3f rx; //roll
//    rx << 1.0f, 0.0f, 0.0f, 0.0f, costheta, -sintheta, 0.0f, sintheta, costheta;
//
//    Eigen::Matrix3f ry; //pitch
//    ry << cosphi, 0.0f, sinphi, 0.0f, 1.0f, 0.0f, -sinphi, 0.0f, cosphi;
//
//    Eigen::Matrix3f rz; //yaw
//    rz << cospsi, -sinpsi, 0.0f, sinpsi, cospsi, 0.0f, 0.0f, 0.0f, 1.0f;
//
    rotationMatrix[0] = vector3(cospsi*cosphi, cospsi*sinphi_sintheta-sinpsi*costheta, cospsi*sinphi_costheta + sinpsi*sintheta);
    rotationMatrix[1] = vector3(sinpsi*cosphi, sinpsi*sinphi_sintheta+cospsi*costheta, sinpsi*sinphi_costheta - cospsi*sintheta);
    rotationMatrix[2] = vector3(-sinphi, cosphi*sintheta, cosphi*costheta);
}

void Engine::rotateTarget(std::vector<vector3> & rotationMatrix, std::vector<vector3> & target, std::vector<vector3> & posed){

    const vector3 * const row0 = &rotationMatrix[0];
    const vector3 * const row1 = &rotationMatrix[1];
    const vector3 * const row2 = &rotationMatrix[2];

    unsigned int count = 0;
    for(auto & vec : target ){
        posed[count] = vector3((*row0).dot(vec),(*row1).dot(vec), (*row2).dot(vec));
        count++;
    }
}


void Engine::translateTarget(vector3 & translation_vector, std::vector<vector3> & target){

    for(auto & vec : target ){
        vec += translation_vector;
    }
}


float Engine::getChiFromFile(std::string filename){

    std::ifstream fs;
    fs.open(filename);
    if(fs.is_open()) {
        fs.seekg(-1,std::ios_base::end);                // go to one spot before the EOF

        if(fs.peek() == '\n') {
            fs.seekg(-1, std::ios_base::cur);
            unsigned int i = fs.tellg();
            for(i;i > 0; i--)
            {
                if(fs.peek() == '\n')
                {
                    //Found
                    fs.get();
                    break;
                }
                fs.seekg(i, std::ios_base::beg);
            }
        }

        std::string lastLine;
        getline(fs,lastLine);                      // Read the current line
        fs.close();

        std::vector<std::string> contents;
        boost::split(contents, lastLine, boost::is_any_of("\t  "), boost::token_compress_on);

        float chi = std::strtof(contents[contents.size()-1].c_str(), nullptr);

        return (std::isnan(chi) ? 1000.0f : chi);
    }

    return 1000.0f;
}


float Engine::closestDistance(std::vector<vector3> & pose){

    float temp, dis=FLT_MAX;

    for (auto & ref_vec : centered_ref_coordinates){

        for(auto & tar_vec : pose){
            temp = (ref_vec - tar_vec).sqlength();
            if (temp < dis){
                dis = temp;
            }
        }
    }
    return std::sqrt(dis);
}

void Engine::writeMergedCoordinates(std::string name, PDBModel & targetModel, std::vector<vector3> & posed){

    std::string residue_index;

    name = name + ".pdb";
    //const char * outputFileName = name.c_str() ;
    //const char * originalPDBFilename = this->filename.c_str();
    FILE * pFile;
    pFile = fopen(name.c_str(), "w");
    vector3 * vec;

    unsigned int totalCoords = reference_model.getTotalCoordinates();

    auto resids = reference_model.getResid();
    std::string chain = "A";

    for(unsigned int n = 0; n<totalCoords; n++){
        residue_index = boost::lexical_cast<std::string>(resids[n]);
        vec = &centered_ref_coordinates[n];
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", n+1, reference_model.getAtomTypeByIndex(n).c_str(), reference_model.getResidueAt(n).c_str(), chain.c_str(), residue_index.c_str(), vec->x, vec->y, vec->z);

    }
    fprintf(pFile,"TER\n");

    /*
     * do target
     */
    totalCoords = targetModel.getTotalCoordinates();
    resids = targetModel.getResid();
    chain = "B";

    for(unsigned int n = 0; n<totalCoords; n++){
        residue_index = boost::lexical_cast<std::string>(resids[n]);
        vec = &posed[n];
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", n+1, targetModel.getAtomTypeByIndex(n).c_str(), targetModel.getResidueAt(n).c_str(), chain.c_str(), residue_index.c_str(), vec->x, vec->y, vec->z);
    }

    fprintf(pFile,"END\n");
    fclose(pFile);

}