#pragma once

#include <iostream>
#include <algorithm>
#include <tuple>
#include <limits>
#include <cfloat>
#include "Definitions.h"

enum class touristAttractionType {none, sight, route};
enum class sightCategory {none, depot, museum, agora};

struct Position {
	double lat, lon;

	Position() {
		lat = 0;
		lon = 0;
	}
	Position(double pLat, double pLon) {
		lat = pLat;
		lon = pLon;
	}
	~Position() {}
};

struct Point {
	int id;
	Position pos;

	Point() {
		id = DEFAULT_POINT_ID;
	}

	Point(int pId, double pLat, double pLon) {
		id = pId;
		pos = Position(pLat, pLon);
	}
	~Point() {}
};

struct TimeWindow {
	double openTime;
	double closeTime;
};

struct BucketActivity {
	int inBucket = 0;
	int inWalk = 0;
	double duration = 0;
};

typedef struct TouristAttraction {
	const std::string id;
	Point point;
	TimeWindow timeWindow;
	const touristAttractionType type;
	const sightCategory category;
	const double visitDuration = 0;
	const int profit;
	double arrTime, waitDuration, startOfVisitTime, depTime, shift, maxShift;
	int route, cluster;
	double minDist;  // default infinite dist to nearest cluster
	int arrPointId, depPointId;
	TouristAttraction* next;
	TouristAttraction* prev;
	BucketActivity bucketActivities[2];


	//default constructor
	TouristAttraction()
		: id(DEFAULT_TA_ID),
		point(Point()),
		visitDuration(.0),
		profit(0),
		timeWindow(TimeWindow{ .0, .0 }),
		arrTime(.0),
		waitDuration(.0),
		startOfVisitTime(.0),
		depTime(.0),
		shift(.0),
		maxShift(.0),
		route(-1),
		cluster(-1),
		minDist(DBL_MAX),
		arrPointId(DEFAULT_POINT_ID),
		depPointId(DEFAULT_POINT_ID),
		type(touristAttractionType::sight),
		category(sightCategory::none),
		next(nullptr), prev(nullptr)
	{
	}
	
	TouristAttraction(Point p)
		: id(DEFAULT_DEPOT_ID),
		point(p),
		visitDuration(.0),
		profit(0),
		timeWindow(TimeWindow{ .0, .0 }),
		arrTime(.0),
		waitDuration(.0),
		startOfVisitTime(.0),
		depTime(.0),
		shift(.0),
		maxShift(.0),
		route(-1),
		cluster(-1),
		minDist(DBL_MAX),
		arrPointId(p.id),
		depPointId(p.id),
		type(touristAttractionType::sight),
		category(sightCategory::none),
		next(nullptr), prev(nullptr)
	{
	}

	TouristAttraction(std::string pId, Point p, double pVisitDuration, int pProfit, double pOpenTime, double pCloseTime)
		: id(pId),
		point(p),
		visitDuration(pVisitDuration),
		profit(pProfit),
		timeWindow(TimeWindow{ pOpenTime, pCloseTime }),
		arrTime(.0),
		waitDuration(.0),
		startOfVisitTime(.0),
		depTime(.0),
		shift(.0),
		maxShift(.0),
		route(-1),
		cluster(-1),
		minDist(DBL_MAX),
		arrPointId(p.id),
		depPointId(p.id),
		type(touristAttractionType::sight),
		category(sightCategory::none),
		next(nullptr), prev(nullptr)
	{
	}

	double euclidean_distance(TouristAttraction ta) { 
		return (ta.point.pos.lat - point.pos.lat) * (ta.point.pos.lat - point.pos.lat) + (ta.point.pos.lon - point.pos.lon) * (ta.point.pos.lon - point.pos.lon);
	}

	virtual void print() {
		std::cout << "default ta print" << std::endl;
	}

	virtual void printMin(){
		std::cout << "default ta printMin" << std::endl;
	}

	virtual TouristAttraction* clone() const {
		std::cout << "default ta clone" << std::endl;
		return new TouristAttraction(*this);
	};

	void reset() {
		arrTime = .0;
		waitDuration = .0;
		startOfVisitTime = .0;
		depTime = .0;
		shift = .0;
		maxShift = .0;
	}
	

	virtual ~TouristAttraction() {} //virtual destructor to ensure our subclasses are correctly dealocated
} TA;

struct Sight : TA {
	using TouristAttraction::TouristAttraction;

	void print() override {
		std::cout << "id: " << id
			<< "\n point.id: " << point.id
			<< "\n visitDuration: " << visitDuration
			<< "\n profit: " << profit
			<< "\n openTime: " << timeWindow.openTime
			<< "\n closeTime: " << timeWindow.closeTime
			<< "\n arrTime: " << arrTime
			<< "\n waitDuration: " << waitDuration
			<< "\n startOfVisitTime: " << startOfVisitTime
			<< "\n depTime: " << depTime
			<< "\n shift: " << shift
			<< "\n maxShift: " << maxShift
			<< "\n route: " << route
			<< "\n cluster: " << cluster
			<< "\n minDist: " << minDist
			<< std::endl;
	}

	void printMin() override{
		std::string out = "id: " + id;
		if(next != nullptr) {
			out += "\nnext:" + next->id;
		} else {
			out += "\nnext: NULL";
		}

		if(prev != nullptr) {
			out += "\nprev:" + prev->id;
		} else {
			out += "\nprev: NULL";
		}

		std::cout << out << std::endl;
	}

	TouristAttraction* clone() const override {
		return new Sight(*this);
	}

};

struct Route : TA {

	Point edges[2];

	Route(std::string pId, Point p, Point p1, Point p2, double pVisitDuration, int pProfit, double pOpenTime, double pCloseTime) :
		TouristAttraction(pId, p, pVisitDuration, pProfit, pOpenTime, pCloseTime), edges{p1, p2}
	{}

	void print() override {
		std::cout << "id: " << id
			<< "\n point.id: " << point.id
			<< "\n edge1.id " << edges[0].id
			<< "\n edge2.id " << edges[1].id
			<< "\n visitDuration: " << visitDuration
			<< "\n profit: " << profit
			<< "\n openTime: " << timeWindow.openTime
			<< "\n closeTime: " << timeWindow.closeTime
			<< "\n arrTime: " << arrTime
			<< "\n waitDuration: " << waitDuration
			<< "\n startOfVisitTime: " << startOfVisitTime
			<< "\n depTime: " << depTime
			<< "\n shift: " << shift
			<< "\n maxShift: " << maxShift
			<< "\n route: " << route
			<< "\n cluster: " << cluster
			<< "\n minDist: " << minDist
			<< std::endl;
	}

	void printMin() override{
		std::string out = "id: " + id;
		if(next != nullptr) {
			out += "\nnext:" + next->id;
		} else {
			out += "\nnext: NULL";
		}

		if(prev != nullptr) {
			out += "\nprev:" + prev->id;
		} else {
			out += "\nprev: NULL";
		}

		std::cout << out << std::endl;
	}

	TouristAttraction* clone() const override {
		return new Route(*this);
	}

};

struct Cluster {
	
	int id, openTime, closeTime;
	Point centroid;
	std::vector<TA*> nodes;

	Cluster() : id(DEFAULT_CLUSTER_ID), openTime(-1), closeTime(-1), nodes(){}

	Cluster(int openTime, int closeTime, std::vector<TA*> attractions, double** ttMatrix) :
		id(DEFAULT_CLUSTER_ID), openTime(openTime), closeTime(closeTime), nodes(attractions) {}

	void print() {
		std::cout << 
			"Cluster: " << id << 
			" - openTme: " << openTime <<
			" - closeTime: " << closeTime << std::endl << "Sights: " << std::endl;
		for (auto p = nodes.begin(); p != nodes.end(); ++p) {
			(*p)->print();
		}
	}

	void setCentroid(int id) {
		double totalLat = .0;
		double totalLon = .0;
		size_t nodesSize = nodes.size();
		for (auto& node : nodes) {
			totalLat += node->point.pos.lat;
			totalLon += node->point.pos.lon;
		}
		centroid = Point(id, totalLat / nodesSize, totalLon / nodesSize);
	}
};

struct Candidate {
	int id;
	int bestPos;
	int bestArrPointId;
	int bestDepPointId;
};


template<class T>
class List {
private:
	void disconnect(T* curr){
		if(curr->next != nullptr){
			curr->next->prev = curr->prev;
		}
		if(curr->prev != nullptr){
			curr->prev->next = curr->next;
		}
		curr->next = nullptr;
		curr->prev = nullptr;
	}

protected:
	T* head;
	T* tail;
	int length;

public:

	List<T>() {
		head = nullptr;
		tail = nullptr;
		length = 0;
	}

	List<T>(T* start) {
		head = start;
		T* curr = head;
		while (curr->next != nullptr) {
			curr = curr->next;
		}
		tail = curr;
	}

	List<T>(T* start, T* last) {
		head = start;
		tail = last;
	}

	T* first() const {
		return head;
	}

	T* second() const {
		if (head->next != nullptr) {
			return head->next;
		}
	}

	T* last() const {
		return tail;
	}

	void print() {
		char type = '_';
		Sight* s = nullptr;
		Route* r = nullptr;
		if (head != nullptr)
		{
			T* curr = head;
			while (curr != nullptr)
			{
				if ((s = dynamic_cast<Sight*>(curr))) {
					type = 's';
				}
				else if ((r = dynamic_cast<Route*>(curr))) {
					type = 'r';
				}
				if (curr == head) {
					std::cout << curr->id << "(" << curr << ")" << "[" << type << "]" << "(head)" << "\t <=> \t";
				}
				else if (curr == tail) {
					std::cout << curr->id << "(" << curr << ")" << "[" << type << "]" << "(tail)";
				}
				else {

					std::cout << curr->id << "(" << curr << ")" << "[" << type << "]" << "\t<=>\t";
				}
				curr = curr->next;
			}

			// std::cout << "<=\ttail(" << tail << ")";
			std::cout << std::endl;
		}
		else
		{
			std::cout << "List is empty." << std::endl;
		}
	}

	void print(std::string msg) {
		std::cout << msg << std::endl;
		char type = '_';
		Sight* s = nullptr;
		Route* r = nullptr;
		if (head != nullptr)
		{
			T* curr = head;
			while (curr != nullptr)
			{
				if ((s = dynamic_cast<Sight*>(curr))) {
					type = 's';
				}
				else if ((r = dynamic_cast<Route*>(curr))) {
					type = 'r';
				}
				if(curr == head){
					std::cout << curr->id << "(" << curr << ")" << "[" << type <<  "]" << "(head)" << "\t <=> \t";
				} else if (curr == tail){
					std::cout << curr->id << "(" << curr << ")" << "[" << type << "]" << "(tail)";
				} else {
					
					std::cout << curr->id << "(" << curr << ")" << "[" << type << "]" << "\t<=>\t";
				}
				curr = curr->next;
			}
			
			// std::cout << "<=\ttail(" << tail << ")";
			std::cout << std::endl;
		}
		else
		{
			std::cout << "List is empty." << std::endl;
		}
	}

	

	T* grabUsingPos(int index) {
		T* curr = head;
		T* temp;
		int pos = 0;
		while (curr != nullptr) {
			if (pos == index) {
				temp = curr;
				curr = curr->next;
				if (temp == head) {
					head = temp->next;
				}
				if (temp == tail) {
					tail = temp->prev;
				}
				if (temp->prev != nullptr) {
					temp->prev->next = temp->next;
				}
				if (temp->next != nullptr) {
					temp->next->prev = temp->prev;
				}
				return temp;
			}
			else {
				pos++;
				curr = curr->next;
			}
		}
		return nullptr;
	}

	//i am not sure if I should do that. It saves some times by not searching for an id, but what if the n isn't inside the current list?
	T* grab(T* n) {
		if (n == head) {
			head = n->next;
		}

		if (n == tail) {
			tail = tail->next;
		}

		if (n->prev != nullptr) {
			n->prev->next = n->next;
		}

		if (n->next != nullptr) {
			n->next->prev = n->prev;
		}
		return n;
	}

	T* grabUsingID(std::string id) {
		T* curr = head;
		T* temp;
		while (curr != nullptr) {
			if (curr->id == id) {
				temp = curr;
				curr = curr->next;
				if (temp == head) {
					head = temp->next;
				}
				if (temp == tail) {
					tail = temp->prev;
				}
				if (temp->prev != nullptr) {
					temp->prev->next = temp->next;
				}
				if (temp->next != nullptr) {
					temp->next->prev = temp->prev;
				}
				return temp;
			}
			else {
				curr = curr->next;
			}
			
		}
		return nullptr;
	}

	void removeFirst(){
		T* curr = head;
		head = curr->next;
		disconnect(curr);
		delete(curr);
	}

	void removeLast(){
		T* curr = tail;
		tail = curr->prev;
		disconnect(curr);
		delete(curr);
	}

	//maybe i could improve this function
	void remove(int index) {
		T* curr = head;
		int pos = 0;
		while (curr != nullptr) {
			if (pos == index) {
				if(curr == head){
					head = curr->next;
				}

				if(curr == tail){
					tail = curr->prev;
				}

				disconnect(curr);
				delete(curr);
			}
			else {
				pos++;
				curr = curr->next;
			}
		}
	}

	void push(T node) {
		T* n = &node;
		if (head == nullptr) {
			// The list is empty
			head = n; 
			head->next = nullptr;
			head->prev = nullptr;
			tail = head;
		}
		else {
			// The list isn't empty
			if (tail == head) {
				// The list has one element
				tail = n;
				tail->next = nullptr;
				tail->prev = head;
				head->next = tail;
			}
			else {
				n->next = nullptr;
				n->prev = tail;
				tail->next = n;
				tail = n;
			}
		}
	}

	void iterate(void (*f)(TA*)) {
		if (head == nullptr) {
			return;
		}
		else {
			T* curr = head;
			while (curr != nullptr) {
				(*f)(curr);
				curr = curr->next;
			}
		}
	
	}
	
	void push(T* ptr) {
		if (head == nullptr) {
			// The list is empty
			ptr->next = nullptr;
			ptr->prev = nullptr;
			head = ptr;
			tail = ptr;
		}
		else {
			ptr->next = nullptr;
			ptr->prev = tail;
			tail->next = ptr;
			tail = ptr;
		}
	}

	void pushNew(T* ref) {
		T* n = ref->clone();
		//*n = *ref;
		if (head == nullptr) {
			// The list is empty
			head = n;
			head->next = nullptr;
			head->prev = nullptr;
			tail = head;
		} else {
			n->next = nullptr;
			n->prev = tail;
			tail->next = n;
			tail = n;
		}
	}

	bool isEmpty() {
		return head == nullptr;
	}

	///* Function to delete the entire linked list */
	//void del()
	//{
	//	/* deref head_ref to get the real head */
	//	Node<T>* curr = head;
	//	Node<T>* next;

	//	while (curr != nullptr)
	//	{
	//		next = current->next;
	//		free(curr);
	//		curr = next;
	//	}

	//	/* deref head_ref to affect the real head back
	//		in the caller. */
	//	*head = NULL;
	//}



	T* get(int index) {
		T* curr = this->head;
		for (int i = 0; i < index; ++i) {
			curr = curr->next;
		}
		return curr;
	}

	T* operator[](int index) {
		return get(index);
	}

	int getLength() {
		int length = 0;
		if (head != NULL)
		{
			T* curr = head;
			while (curr != nullptr)
			{
				length++;
				curr = curr->next;
			}
		}

		return length;
	}

	List<T>* splitAtIndex(int index) {
		if (index < 0) {
			index = getLength() - index;
		}
		T* curr = head;
		std::tuple<List<T>, List<T>> parts;
		int pos = 0;
		while (curr != nullptr) {
			if (pos >= index) {
				std::get<1>(parts).pushNew(curr);
			}
			else {
				std::get<0>(parts).pushNew(curr);
			}
			
			pos++;
			curr = curr->next;
		}

		return parts;
	}

	void setLength() {
		this->length = this->getLength();
	}

	void append(List<T> list)
	{
		if (head == nullptr) { // list is empty
			head = list.head;
			tail = list.tail;
		}
		else {
			tail->next = list.head;
			list.head->prev = tail;
			tail = list.tail;
		}
		list.head = nullptr;
		list.tail = nullptr;
	}

	void empty() {
		T* curr = head;
		T* del;
		while (curr != nullptr) {
			del = curr;
			curr = curr->next;
			delete(del);
		}

		head = nullptr;
		tail = nullptr;
	}
};

typedef class TouristAttractionList :  public List<TA> {

public:

	TouristAttractionList() {

	}

	TouristAttractionList(TA* start) {
		head = start;
		TA* curr = head;

		while (curr->next != nullptr) {
			curr = curr->next;
		}
		tail = curr;
	}

	TouristAttractionList(std::vector<TA*> v) {
		for (auto& i : v) {
			this->push(i);
		}
	}

	void printMin(std::string msg) {
		std::cout << msg << std::endl;
		TA* curr = head;
		while (curr != nullptr) {
			std::cout << curr->id << "\t";
			curr = curr->next;
		}
		std::cout << std::endl;
	}

	//insert ta n right before index
	void insertAt(TA* n, int index, int arrPointId, int depPointId, std::vector<std::vector<double>> ttMatrix) {

		if (head == nullptr) {  // The list is empty
			head = new TA;
			head = n;
			tail = head;
		}
		else { // The list isn't empty
			if (tail == head) { // The list has one element
				head->next = n;
				n->prev = head;
				tail = n;
			}
			else {
				// The list has more than one element
				TA* curr = head;
				TA* temp;
				int pos = 0;

				while (curr != nullptr) {
					if (pos == index) {
						if (curr == head) {
							std::exit(1);
						} else {
							temp = curr->prev;
							n->shift = ttMatrix[temp->depPointId][arrPointId] + n->waitDuration + n->visitDuration + ttMatrix[depPointId][curr->arrPointId] - ttMatrix[temp->depPointId][curr->arrPointId];
							n->arrPointId = arrPointId;
							n->depPointId = depPointId;
							
						}

						if (curr->prev != nullptr) {
							curr->prev->next = n;
						}
						n->prev = curr->prev;
						curr->prev = n;
						n->next = curr;
						return;
					}
					else {
						curr = curr->next;
						pos++;
					}
				}
			}
		}
	}

	std::vector<TA> toVec() {
		std::vector<TA> v;
		TA* curr = head;
		while (curr != nullptr) {
			v.push_back(*curr);
			curr = curr->next;
		}
		return v;
	}

	TouristAttractionList copy() {
		TouristAttractionList copy;
		TA* curr = head;
		while (curr != nullptr) {
			copy.pushNew(curr);
			curr = curr->next;
		}
		return copy;
	}

	TouristAttractionList copyPart(int startIndex, int endIndex) {

		TouristAttractionList copy;
		TA* curr = head;
		bool reachedEndIndex = false;
		int pos = 0;

		if (startIndex < 0 || endIndex < 0) {
			size_t length = getLength();
			if (startIndex < 0) {
				startIndex += length;
			}

			if (endIndex < 0) {
				endIndex += length;
			}
		}

		if (startIndex > endIndex) {
			std::cerr << "Invalid startIndex or endIndex" << std::endl;
			std::exit(1);
		}

		while (curr != nullptr) {
			if (pos >= startIndex) {
				copy.pushNew(curr);

				if (pos == endIndex) {
					reachedEndIndex = true;
					break;
				}
			}

			pos++;
			curr = curr->next;
		}

		if (!reachedEndIndex) {
			std::cout << "Didn't found the end index so returned copied until the end of the list" << std::endl;
		}

		return copy;
	}

	//removePart removes a part of the list
	// si = start index
	// ei = end index 
	TouristAttractionList removePart(int si, int ei)
	{

		if(si > ei) {
			std::cerr << "Invalid arguments" << std::endl;
			exit(1);
		}

		TA* curr = head;
		TA* temp;
		int pos = 0;

		while (curr != nullptr) {
			if (pos == si) {
				temp = curr; //temp points to the first node that gets deleted
				break;
			}
			else {
				pos++;
				curr = curr->next;
			}

		}

		while(curr != nullptr){
			if(pos == ei){

				// if (curr == tail) {
				// 	curr = tail->prev;
				// }
				std::cout << "Will remove part from " << temp->id << " to " << curr->id << " ( " << si << " " << ei << " of total route length " <<  this->getLength() << " ) " << std::endl;

				// if (head == temp) {
				// 	head = curr->next;
				// }

				if (temp->prev != nullptr) {
					temp->prev->next = curr->next;
				}

				if (curr->next != nullptr) {
					curr->next->prev = temp->prev;
				}

				temp->prev = nullptr;
				curr->next = nullptr;


				//TODO: delete nodes
				TouristAttractionList removed(temp);
				return removed;

			} else {
				pos++;
				curr = curr->next;
			}
		}



		return *this;
	}

	//grabPart grabs a part of the list
	// si = start index
	// ei = end index 
	TouristAttractionList grabPart(int si, int ei)
	{

		if(si > ei) {
			std::cerr << "Invalid arguments" << std::endl;
			exit(1);
		}

		TA* curr = head;
		TA* temp;
		int pos = 0;

		while (curr != nullptr) {
			if (pos == si) {
				temp = curr; //temp points to the first node that gets deleted
				break;
			}
			else {
				pos++;
				curr = curr->next;
			}

		}

		while(curr != nullptr){
			if(pos == ei){

				if (temp->prev != nullptr) {
					temp->prev->next = curr->next;
				}

				if (curr->next != nullptr) {
					curr->next->prev = temp->prev;
				}

				temp->prev = nullptr;
				curr->next = nullptr;

				TouristAttractionList part(temp);
				return part;

			} else {
				pos++;
				curr = curr->next;
			}
		}
		return *this;
	}

	void removeNodesWithId(std::string id) {
		TA* curr = head;
		TA* temp;
		while (curr != nullptr) {
			if (curr->id == id) {
				temp = curr;
				curr = curr->next;
				if (temp == head) {
					head = temp->next;
				}
				if (temp == tail) {
					tail = temp->prev;
				}
				if (temp->next != nullptr) {
					temp->next->prev = temp->prev;
				}
				if (temp->prev != nullptr) {
					temp->prev->next = temp->next;
				}
				temp->next = nullptr;
				temp->prev = nullptr;
				delete(temp);
			}
			else {
				curr = curr->next;
			}			
		}
	}

	std::tuple<TouristAttractionList, TouristAttractionList> splitAtIndex(int index) {
		TA* curr = head;
		std::tuple<TouristAttractionList, TouristAttractionList> lists;
		int pos = 0;

		if (index < 0) {
			index = getLength() - index;
		}

		while (curr != nullptr) {
			if (pos >= index) {
				std::get<1>(lists).pushNew(curr);
			}
			else {
				std::get<0>(lists).pushNew(curr);
			}

			pos++;
			curr = curr->next;
		}

		return lists;
	}

	std::tuple<TouristAttractionList, TouristAttractionList> splitOnDepTime(double depTime) {
		TA* curr = head;
		std::tuple<TouristAttractionList, TouristAttractionList> lists;
		while (curr != nullptr) {
			if (curr->depTime <= depTime) {
				std::get<0>(lists).pushNew(curr);
			}
			else {
				std::get<1>(lists).pushNew(curr);
			}
			curr = curr->next;
		}
		return lists;
	}

	double getDuration() {
		if (head == nullptr) {
			return 0;
		}
		return tail->depTime - head->arrTime;
	}

	int collectProfit() {
		TA* curr = head;
		int totalProfit = 0;
		while (curr != nullptr) {
			totalProfit += curr->profit;
			curr = curr->next;
		}
		return totalProfit;
	}

	bool isValid(bool check_window, double Tmax, std::vector<std::vector<double>>& ttMatrix ) {
		TA* curr = head;
		if (curr == nullptr) {
			std::cout << "List is empty" << std::endl;
			return false;
		}
		while (curr != nullptr) {
			if (curr->startOfVisitTime != curr->arrTime + curr->waitDuration) {
				std::cout << "[" << curr->id << "] StartOfVisitTime(" << curr->startOfVisitTime << ") != arrTime(" << curr->arrTime << ")" << std::endl;
				return false;
			}

			if (curr->depTime != curr->startOfVisitTime + curr->visitDuration) {
				std::cout << "[" << curr->id << "] depTime(" << curr->depTime << ") != startOfVisitTime(" << curr->startOfVisitTime << ") + visitDuration(" << curr->visitDuration << ")" << std::endl;
				return false;
			}

			if (check_window && curr->depTime > curr->timeWindow.closeTime) {
				std::cout << "[" << curr->id << "] depTime(" << curr->depTime << ") > closeTime(" << curr->timeWindow.closeTime << ")" << std::endl;
				return false;
			}

			if (curr->next != nullptr) {
				if (curr->next->arrTime != curr->depTime + ttMatrix[curr->depPointId][curr->next->arrPointId]) {
					std::cout << "[" << curr->id << "] next->arrTime(" << curr->next->arrTime << ") != depTime(" << curr->depTime << ") + timeTravel[" << curr->depPointId << "][" << curr->next->arrPointId << "](" << ttMatrix[curr->depPointId][curr->next->arrPointId] << ")" << std::endl;
					return false;
				}
			}

			curr = curr->next;
		}

		return true;
	}
} ListTA; 

using Walk = ListTA;

