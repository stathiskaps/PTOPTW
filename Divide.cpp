#include "Divide.h"

Divide::Divide(){

}

Divide::Divide(std::vector<TA*> attractions){
    mAttractions = attractions;
}

Divide::~Divide(){

}

int Divide::calcMeanTimePoint(){
    int sum = 0;
    double meanVal = .0;
    double dist;
    double min = DBL_MAX;
    int timePoint = -1;
    for(auto &a: mAttractions){
        sum += a->timeWindow.openTime;
        sum += a->timeWindow.closeTime;
    }
    
    meanVal = sum / (mAttractions.size() * 2);

    std::cout << "meanVal:" << meanVal << std::endl;



    for (auto &a: mAttractions){
        dist = abs(a->timeWindow.openTime - meanVal);
        if(dist < min){
            min = dist;
            timePoint = a->timeWindow.openTime;
            std::cout << "foundNewMin at " << timePoint << std::endl;
        }
        dist = abs(a->timeWindow.closeTime - meanVal);
        if(dist < min){
            min = dist;
            timePoint = a->timeWindow.closeTime;
            std::cout << "foundNewMin at " << timePoint << std::endl;
        }
    }

    return timePoint;
}

bool Divide::compareByProfit(const TA* a, const TA* b) {
    return a->profit < b->profit;
}

std::vector<OP> Divide::exec(){
    int timePoint = calcMeanTimePoint();
    std::cout << "timePoint: " << timePoint << std::endl;
    std::vector<OP> example;

    std::vector<std::string> A, B;



    std::vector<TA*> attractions = mAttractions;
    std::sort(attractions.begin(), attractions.end(), compareByProfit);

    for(auto &ta: attractions){
        std::cout << ta->profit << " " << std::endl;
    }

    for(auto &ta : attractions){
        int leftDuration = timePoint - ta->timeWindow.openTime;
        int rightDuration = ta->timeWindow.closeTime - timePoint;

        if(leftDuration > rightDuration){
            A.push_back(ta->id);
        } else {
            std::cout <<"come oonnn" << std::endl;
            B.push_back(ta->id);
        }
    }

    std::cout << "Printing A list" << std::endl;
    std::cout << A.size() << std::endl;
    for(auto &id: A){
        std::cout << id << " ";
    }

    std::cout << std::endl;

    std::cout << "Printing B list" << std::endl;
    std::cout << B.size() << std::endl;
    for(auto &id: B){
        std::cout << id << " ";
    }

    return example;

}