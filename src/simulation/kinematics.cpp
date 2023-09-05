#include "simulation/kinematics.h"

#include "Eigen/Dense"
#include <iostream>
#include "acclaim/bone.h"
#include "util/helper.h"
#include <algorithm>
#include <stack>
#include <map>
#include <cmath>

namespace kinematics {

void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
    // TODO (FK)
    // Same as HW2
    // Hint:
    //   1. If you don't use `axis` in this function, you can copy-paste your code

    bone->start_position = Eigen::Vector4d::Zero();
    bone->end_position = Eigen::Vector4d::Zero();
    bone->rotation = Eigen::Matrix4d::Zero();
    bone->start_position = Eigen::Vector4d::Zero();
    bone->end_position = Eigen::Vector4d::Zero();
    bone->rotation = Eigen::Matrix4d::Zero();

    std::map<int, bool> visit;
    std::stack<acclaim::Bone*> q;
    visit[bone->idx] = true;
    bone->start_position = posture.bone_translations[bone->idx];
    bone->end_position = posture.bone_translations[bone->idx];
    q.push(bone);

    while (!q.empty()) {
        acclaim::Bone* t = q.top();
        q.pop();
        if (t->idx != 0) {
            t->start_position = t->parent->end_position;
            Eigen::Quaterniond localRotation = util::rotateDegreeZYX(posture.bone_rotations[t->idx].x(), posture.bone_rotations[t->idx].y(), posture.bone_rotations[t->idx].z());
            Eigen::Affine3d rot = t->rot_parent_current * Eigen::Affine3d(localRotation);
            for (acclaim::Bone* itr = t->parent; itr != nullptr; itr = itr->parent) {
                Eigen::Quaterniond parentLocal = util::rotateDegreeZYX(posture.bone_rotations[itr->idx].x(), posture.bone_rotations[itr->idx].y(), posture.bone_rotations[itr->idx].z());
                rot = itr->rot_parent_current * parentLocal * rot;
            }

            Eigen::Vector4d dir = t->dir * t->length;
            t->end_position = rot.matrix() * dir + posture.bone_translations[t->idx] + t->start_position;

            t->rotation = rot;
        }
        for (acclaim::Bone* itr = t->child; itr != nullptr; itr = itr->sibling) {
            if (visit.find(itr->idx) == visit.end()) {
                visit[itr->idx] = true;
                q.push(itr);
            }
        }
    }
}

Eigen::VectorXd pseudoInverseLinearSolver(const Eigen::Matrix4Xd& Jacobian, const Eigen::Vector4d& target) {
    // TODO (find x which min(| jacobian * x - target |))
    // Hint:
    //   1. Linear algebra - least squares solution
    //   2. https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Construction
    // Note:
    //   1. SVD or other pseudo-inverse method is useful
    //   2. Some of them have some limitation, if you use that method you should check it.
    Eigen::VectorXd deltatheta(Jacobian.cols());
    deltatheta.setZero();
    // J(theta) = target
    // theta = (J^+)*V
    // SVD : J = U E V^T, J^+ = V E^+ U^T 
    if (Jacobian.isZero(0))
        return deltatheta;
    //Eigen::EigenSolver<Eigen::MatrixXd> eigens(Jacobian.transpose() * Jacobian);
    //Eigen::MatrixXcd x = eigens.eigenvectors();
    //Eigen::MatrixXd v(Jacobian.cols(), Jacobian.cols());
    //v.setZero();
    //for (int i = 0; i < x.cols(); i++)
    //    v.col(i) = x.col(i).real().normalized();

    //Eigen::MatrixXd sigma(Jacobian.rows(), Jacobian.cols());
    //for (int i = 0; i < eigens.eigenvalues().cols(); i++)
    //    sigma.col(i)[i] = sqrt(eigens.eigenvalues()[i].real());
    //Eigen::MatrixXd u(Jacobian.rows(), Jacobian.rows());
    //for (int i = 0; i < v.cols(); i++)
    //    sigma.col(i) = Jacobian * v.col(i) / sigma.col(i)[i];

    //// J = u sigma v^T, J^+ = V E^+ U^T
    //Eigen::Matrix4Xd IJacobian = v * sigma.inverse() * u.transpose();
    //deltatheta = IJacobian * target;
    Eigen::MatrixXd pinv =  Jacobian.completeOrthogonalDecomposition().pseudoInverse();
    deltatheta = pinv *  target;
    
    return deltatheta;
}

/**
 * @brief Perform inverse kinematics (IK)
 *
 * @param target_pos The position where `end_bone` will move to.
 * @param start_bone This bone is the last bone you can move while doing IK
 * @param end_bone This bone will try to reach `target_pos`
 * @param posture The original AMC motion's reference, you need to modify this
 *
 * @return True if IK is stable (HW3 bonus)
 */
bool inverseJacobianIKSolver(const Eigen::Vector4d& target_pos, acclaim::Bone* start_bone, acclaim::Bone* end_bone,
    acclaim::Posture& posture) {
    constexpr int max_iteration = 1000;
    constexpr double epsilon = 1E-3;
    constexpr double step = 0.1;
    // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is the root.
    acclaim::Bone* root_bone = start_bone - start_bone->idx;
    // TODO
    // Perform inverse kinematics (IK)
    // HINTs will tell you what should do in that area.
    // Of course you can ignore it (Any code below this line) and write your own code.
    acclaim::Posture original_posture(posture);

    size_t bone_num = 0;
    std::vector<acclaim::Bone*> boneList;

    // TODO
    // Calculate number of bones need to move to perform IK, store in `bone_num` 
    // a.k.a. how may bones from end_bone to its parent then to start_bone (include both start_bone and end_bone)
    // Store the bones need to move to perform IK into boneList
    // Hint:
    //   1. Traverse from end_bone to start_bone is easier than start to end (since there is only 1 parent)
    //   2. If start bone is not reachable from end. Go to root first.
    // Note:
    //   1. Both start_bone and end_bone should be in the list
    acclaim::Bone* current = end_bone;
    while (current != start_bone && current != root_bone)
    {
        //std::cout << current << std::endl;
        bone_num++;
        boneList.emplace_back(current);
        current = current->parent;
    }
    if (current == root_bone) {
        current = start_bone;
        while (current != root_bone)
        {
            bone_num++;
            boneList.emplace_back(current);
            current = current->parent;
        }
    }
    bone_num++;
    boneList.emplace_back(current);

    Eigen::Matrix4Xd Jacobian(4, 3 * bone_num);
    Jacobian.setZero();
    for (int iter = 0; iter < max_iteration; ++iter) {
        forwardSolver(posture, root_bone);
        Eigen::Vector4d desiredVector = target_pos - end_bone->end_position;
        if (desiredVector.norm() < epsilon) {
            break;
        }
        // TODO (compute jacobian)
        //   1. Compute arm vectors
        //   2. Compute jacobian columns, store in `Jacobian`
        // Hint:
        //   1. You should not put rotation in jacobian if it doesn't have that DoF.
        //   2. jacobian.col(/* some column index */) = /* jacobian column */
        for (long long i = 0; i < bone_num; i++) {
            int startCol = i * 3;
            Eigen::Vector3d tempVector = (end_bone->end_position - boneList[i]->start_position).head(3);
            Eigen::Vector3d head = boneList[i]->rotation.matrix().col(0).head(3);
            Eigen::Vector3d temp = head.cross(tempVector);
            Eigen::Vector4d result = Eigen::Vector4d(temp.x(), temp.y(), temp.z(), 0);
            if (boneList[i]->dofrx)
                Jacobian.col(startCol) = result.normalized();

            head = boneList[i]->rotation.matrix().col(1).head(3);
            temp = head.cross(tempVector);
            result = Eigen::Vector4d(temp.x(), temp.y(), temp.z(), 0);
            if (boneList[i]->dofry)
                Jacobian.col(startCol + 1) = result.normalized();

            head = boneList[i]->rotation.matrix().col(2).head(3);
            temp = head.cross(tempVector);
            result = Eigen::Vector4d(temp.x(), temp.y(), temp.z(), 0);
            if (boneList[i]->dofrz)
                Jacobian.col(startCol + 2) = result.normalized();
        }
        Eigen::VectorXd deltatheta = step * pseudoInverseLinearSolver(Jacobian, desiredVector);

        // TODO (update rotation)
        //   Update `posture.bone_rotation` (in euler angle / degrees) using deltaTheta
        // Hint:
        //   1. You can ignore rotation limit of the bone.
        // Bonus:
        //   1. You cannot ignore rotation limit of the bone.
        for (long long i = 0; i < bone_num; i++) {
            double newValue = posture.bone_rotations[boneList[i]->idx].x() + step * deltatheta(i * 3) * 180 / 3.14;
            if (newValue > boneList[i]->rxmin && newValue < boneList[i]->rxmax)
                posture.bone_rotations[boneList[i]->idx].x() = newValue;
            newValue = posture.bone_rotations[boneList[i]->idx].y() + step * deltatheta(i * 3 + 1) * 180 / 3.14;
            if (newValue > boneList[i]->rymin && newValue < boneList[i]->rymax)
                posture.bone_rotations[boneList[i]->idx].y() = newValue;
            newValue = posture.bone_rotations[boneList[i]->idx].z() + step * deltatheta(i * 3 + 2) * 180 / 3.14;
            if (newValue > boneList[i]->rzmin && newValue < boneList[i]->rzmax)
                posture.bone_rotations[boneList[i]->idx].z() = newValue;
        }
    }
    // TODO (Bonus)
    // Return whether IK is stable (i.e. whether the ball is reachable) and let the skeleton not swing its hand in the air
    if ((end_bone->end_position - target_pos).norm() < epsilon)
        return true;
    else {
        posture = original_posture;
        return false;
    }
}
}  // namespace kinematics
