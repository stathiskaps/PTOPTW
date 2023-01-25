#include "Solution.h"

Solution::Solution(){}
Solution::~Solution(){
    m_unvisited.clear();
    for(auto& walk:m_walks){
        walk.clear();
    }
}
Solution::Solution(List<TA> unvisited): m_unvisited(unvisited) {}

//todo: check what explicit is
Solution::Solution(TA start, TA end, List<TA> unvisited, double startTime, double endTime, int walksNum) : m_unvisited(unvisited) {
    for (int i = 0; i < walksNum; ++i) {
        List<TA> walk{ start, end };
        walk.front().depTime = startTime;
        walk.front().timeWindow = TimeWindow{ startTime, endTime };
        walk.back().timeWindow = TimeWindow{ startTime, endTime };
        walk.back().maxShift = endTime - startTime;
        m_walks.push_back(walk);
    }
}

Solution::Solution(TA start, TA end, List<TA> unvisited, size_t walks) : m_unvisited(unvisited){
    for(size_t i = 0; i < walks; ++i){
         m_walks.push_back(List<TA>{start, end});
    }
}

void Solution::print(std::string tag, const bool verbose){
    std::cout << tag << std::endl;

    m_unvisited.print("Unvisited", false);
    for(size_t i = 0; i < m_walks.size(); ++i) {
        m_walks[i].print("Walk "+std::to_string(i), verbose);
    }
}

int Solution::getScores() const{
    int sum{}; 
    for (auto& w : m_walks) {
        for (auto& p : w) {
            sum += p.profit;
        }
    }
    return sum;
}

int Solution::getVisits(){
    int visits{};
    for (auto& w : m_walks) visits += w.size() - 1;
    return visits;
}

int Solution::getMinWalkSize() const {
    int min_size = INT_MAX;
    for (auto& walk : m_walks) {
        if (walk.size() < min_size) min_size = walk.size();
    }
    return min_size;
}

void Solution::draw(std::string tag, std::string color){
    std::string filename = tag+".svg";
    Bounds bounds;
    for(List<TA>::iterator ta_it = m_unvisited.begin(); ta_it != m_unvisited.end(); ++ta_it){
        if(ta_it.iter->data.point.pos.lat < bounds.minLat){
            bounds.minLat = ta_it.iter->data.point.pos.lat;
        }
        if(ta_it.iter->data.point.pos.lat > bounds.maxLat){
            bounds.maxLat = ta_it.iter->data.point.pos.lat;
        }
        if(ta_it.iter->data.point.pos.lon < bounds.minLon){
            bounds.minLon = ta_it.iter->data.point.pos.lon;
        }
        if(ta_it.iter->data.point.pos.lon > bounds.maxLon){
            bounds.maxLon = ta_it.iter->data.point.pos.lon;
        }
    }

    for(auto& walk : m_walks){
        for(List<TA>::iterator ta_it = walk.begin(); ta_it != walk.end(); ++ta_it){
            if(ta_it.iter->data.point.pos.lat < bounds.minLat){
                bounds.minLat = ta_it.iter->data.point.pos.lat;
            }
            if(ta_it.iter->data.point.pos.lat > bounds.maxLat){
                bounds.maxLat = ta_it.iter->data.point.pos.lat;
            }
            if(ta_it.iter->data.point.pos.lon < bounds.minLon){
                bounds.minLon = ta_it.iter->data.point.pos.lon;
            }
            if(ta_it.iter->data.point.pos.lon > bounds.maxLon){
                bounds.maxLon = ta_it.iter->data.point.pos.lon;
            }
        }
    }


	const double radius = 2;

    // Open an output stream to write the SVG file
    std::ofstream out(filename);



    // Write the SVG file header
    out << "<svg viewBox=\"" << bounds.minLat - radius << " " << bounds.minLon - radius << " " << bounds.maxLat - bounds.minLat + radius*2 << " " << bounds.maxLon - bounds.minLon + radius*2 << "\" xmlns=\"http://www.w3.org/2000/svg\">" << std::endl;

	// Write the routes as lines in the SVG file
    for(auto& walk : m_walks){
        List<TA>::iterator left, right;
        for(left = walk.begin(), right = walk.begin() + 1; right != walk.end(); ++left, ++right){
            out << "<line x1=\"" << left.iter->data.point.pos.lat << "\" y1=\"" << left.iter->data.point.pos.lon <<
             "\" x2=\"" << right.iter->data.point.pos.lat << "\" y2=\"" << right.iter->data.point.pos.lon << "\" style=\"stroke:rgb(66,66,66);stroke-width:0.5\" />" << std::endl;
        }
    }

    //Write nodes
    for(List<TA>::iterator ta_it = m_unvisited.begin(); ta_it != m_unvisited.end(); ++ta_it){
        out << "<circle cx=\"" << ta_it.iter->data.point.pos.lat << "\" cy=\"" << ta_it.iter->data.point.pos.lon  << "\" r=\"" << radius << "\" />" << std::endl;
    }
    for(auto& walk : m_walks){
        for(List<TA>::iterator ta_it = walk.begin(); ta_it != walk.end(); ++ta_it){
            out << "<circle cx=\"" << ta_it.iter->data.point.pos.lat << "\" cy=\"" << ta_it.iter->data.point.pos.lon  << "\" r=\"" << radius << "\" />" << std::endl;
        }
    }

    //Write ids
    for(List<TA>::iterator ta_it = m_unvisited.begin(); ta_it != m_unvisited.end(); ++ta_it){
        out << "<text x=\""<< ta_it.iter->data.point.pos.lat << "\" y=\""<< ta_it.iter->data.point.pos.lon << "\" text-anchor=\"middle\" font-size=\"2px\" fill=\"white\" alignment-baseline=\"middle\">" << ta_it.iter->data.id  << "</text>" << std::endl;
    }
    for(auto& walk : m_walks){
        for(List<TA>::iterator ta_it = walk.begin(); ta_it != walk.end(); ++ta_it){
            out << "<text x=\""<< ta_it.iter->data.point.pos.lat << "\" y=\""<< ta_it.iter->data.point.pos.lon << "\" text-anchor=\"middle\" font-size=\"2px\" fill=\"white\" alignment-baseline=\"middle\">" << ta_it.iter->data.id  << "</text>" << std::endl;
        }
    }

    // Write the SVG file footer
    out << "</svg>" << std::endl;

    // Close the output stream
    out.close();
}


