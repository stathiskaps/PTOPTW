#pragma once

#include <iterator>
#include <initializer_list>
#include <algorithm>
#include <iostream>
#include <type_traits>
#include <vector>
#include <string>

// Type trait helper to identify iterators --------------------------------------------------------
template<typename T, typename = void>
struct is_iterator { static constexpr bool value = false; };
template<typename T>
struct is_iterator<T, typename std::enable_if<!std::is_same<typename std::iterator_traits<T>::value_type, void>::value>::type> {
    static constexpr bool value = true;
};

// The List class ---------------------------------------------------------------------------------
template <typename T>
class List {

    // Sub class for a Node -----------
    struct Node {
        T data{};
        Node* next{};
        Node* previous{};
        Node() {}
        Node(Node* const n, Node* const p) : next(n), previous(p) {}
        Node(Node* const n, Node* const p, const T& d) : next(n), previous(p), data(d) {}
    };

    

    // Protected CustomList data and functions --------
    Node* head{};
    size_t numberOfElements{};

protected:

    void init() { head = new Node(); head->next = head; head->previous = head; numberOfElements = 0; }

public:
    struct iterator;    // Forward declaration

    // Constructor --------------------
    List() { init(); }
    explicit List(const size_t count) { init(); insert(begin(), count); }
    explicit List(const size_t count, const T& value) { init(); insert(begin(), count, value); };
    explicit List(const std::vector<T> v) { init(); insert(begin(), v.begin(), v.end()); }
    explicit List(const std::vector<T*> v) { init(); for (auto& i : v) push_back(*i); }
    template <typename Iter>
    List(const Iter& first, const Iter& last) { init(); insert(begin(), first, last); }
    List(const List& other) { init(), insert(begin(), other.begin(), other.end()); };

    List(List&& other) : head(other.head), numberOfElements(other.numberOfElements) { other.init(); }
    List(const std::initializer_list<T>& il) { init(); insert(begin(), il.begin(), il.end()); }
    template <int N> List(T(&other)[N]) { init(); insert(begin(), std::begin(other), std::end(other)); }
    template <int N> List(const T(&other)[N]) { init(); insert(begin(), std::begin(other), std::end(other)); }


    // Assignment ---------------------
    List& operator =(const List& other) { clear(); insert(begin(), other.begin(), other.end()); return *this; }
    List& operator =(List&& other) { clear(); head = other.head; numberOfElements = other.numberOfElements; other.init(); return *this; }
    List& operator =(const std::initializer_list<T>& il) { clear(); insert(begin(), il.begin(), il.end()); return *this; }
    template <int N> List& operator =(const T(&other)[N]) { clear(); insert(begin(), std::begin(other), std::end(other)); return *this; }
    template <int N> List& operator =(T(&other)[N]) { clear(); insert(begin(), std::begin(other), std::end(other)); return *this; }

    template <typename Iter> void assign(const Iter& first, const Iter& last) { clear(); insert(begin(), first, last); }
    template <int N> void assign(const T(&other)[N]) { clear(); insert(begin(), std::begin(other), std::end(other)); }
    template <int N> void assign(T(&other)[N]) { clear(); insert(begin(), std::begin(other), std::end(other)); }
    void assign(const size_t count, const T& value) { clear(); insert(begin(), count, value); }
    void assign(const std::initializer_list<T>& il) { clear(); insert(begin(), il.begin(), il.end()); }

    // Destructor ---------------------
    ~List() { destroy(); }

    // Element Access -----------------
    T& front() { return *begin(); }
    T& back() { return *(--end()); }

    // Iterators ----------------------
    iterator begin() const { return iterator(head->next, head); }
    iterator end() const { return iterator(head, head); }

    // Capacity -----------------------
    inline size_t size() const { return numberOfElements; }
    inline bool empty() const { return size() == 0; }

    // Modifiers ----------------------
    List& disc(iterator& first, iterator& last) {
        if (head == first) {
            head = last->next;
        } 
        first->prev->next = last->next;
        last->next->prev = first->prev;
        first->prev = nullptr;
        last->next = nullptr;
    }

    void print(std::string tag, bool verbose) const{
        std::cout << tag << ": " <<std::endl;
        std::cout << "TA:\t";
        for(iterator it = begin(); it != end(); ++it){
            std::cout << it.iter->data.id << "\t";
        }
        std::cout << std::endl;

        if(!verbose){
            return;
        }

        std::cout << "arrTime:\t";
        for(iterator it = begin(); it != end(); ++it){
            std::cout << it.iter->data.arrTime << "\t";
        }
        std::cout << std::endl;

        std::cout << "openTime\t";
        for(iterator it=begin(); it!=end(); ++it){
            std::cout << it.iter->data.timeWindow.openTime << "\t";
        }
        std::cout << std::endl;

        std::cout << "depTime:\t";
        for(iterator it = begin(); it != end(); ++it){
            std::cout << it.iter->data.depTime << "\t";
        }
        std::cout << std::endl;

        std::cout << "closeTime\t";
        for(iterator it=begin(); it!= end(); ++it){
            std::cout << it.iter->data.timeWindow.closeTime << "\t";
        }
        std::cout << std::endl;

        std::cout << "shift:\t";
        for(iterator it = begin(); it != end(); ++it){
            std::cout << it.iter->data.shift << "\t";
        }
        std::cout << std::endl;
        
        std::cout << "maxShift:\t";
        for(iterator it = begin(); it != end(); ++it){
            std::cout << it.iter->data.maxShift << "\t";
        }
        std::cout << std::endl;

        std::cout << "waitDur:\t";
        for(iterator it = begin(); it != end(); ++it){
            std::cout << it.iter->data.waitDuration << "\t";
        }
        std::cout << std::endl;

        std::cout << "(size: " << size() << ")" << std::endl;
        std::cout << std::endl;
    }

    void clear() {
        for (Node* nextNode{}, * currentNode(head->next); currentNode != head; currentNode = nextNode) {
            nextNode = currentNode->next;
            delete currentNode;
        }
        delete head;
        init();
    }

    std::pair<List<T>, List<T>> split(iterator it) {
        List<T> leftPart, rightPart;

        if(it == end()) {
            leftPart = List<T>(begin(), it);
            rightPart = List<T>();
        } else if(it == begin()) {
            leftPart = List<T>();
            rightPart = List<T>(it, end()-1);
        } else {
            leftPart = List<T>(begin(), it);
            rightPart = List<T>(it, end());
        }

        return {leftPart, rightPart};
    }

    void destroy() {
        for (Node* nextNode{}, * currentNode(head->next); currentNode != head; currentNode = nextNode) {
            nextNode = currentNode->next;
            delete currentNode;
        }
        delete head;
    }

    iterator insert(const iterator& insertBeforePosition, const T& value) {
        Node* nodeInsertBeforePosition = insertBeforePosition.iter;
        Node* newNode = new Node(nodeInsertBeforePosition, nodeInsertBeforePosition->previous, value);
        nodeInsertBeforePosition->previous = newNode;
        (newNode->previous)->next = newNode;
        ++numberOfElements;
        return iterator(newNode, head);
    }
    iterator insert(const iterator& insertBeforePosition) {
        Node* nodeInsertBeforePosition = insertBeforePosition.iter;
        Node* newNode = new Node(nodeInsertBeforePosition, nodeInsertBeforePosition->previous);
        nodeInsertBeforePosition->previous = newNode;
        (newNode->previous)->next = newNode;
        ++numberOfElements;
        return iterator(newNode, head);
    }
    template <class Iter, std::enable_if_t<is_iterator<Iter>::value, bool> = true>
    iterator insert(const iterator& insertBeforePosition, const Iter& first, const Iter& last) {
        iterator result(insertBeforePosition.iter, head);
        if (first != last) {
            result = insert(insertBeforePosition, *first);
            Iter i(first);
            for (++i; i != last; ++i)
                insert(insertBeforePosition, *i);
        }
        return result;
    }
    iterator insert(const iterator& insertBeforePosition, const size_t& count, const T& value) {
        iterator result(insertBeforePosition.iter, head);
        if (count != 0u) {
            result = insert(insertBeforePosition, value);
            for (size_t i{ 1u }; i < count; ++i)
                insert(insertBeforePosition, value);
        }
        return result;
    }
    iterator insert(const iterator& insertBeforePosition, const std::initializer_list<T>& il) {
        return insert(insertBeforePosition, il.begin(), il.end());
    }

    iterator disconnect(const iterator& posToDisconnect) {
        Node* nodeToDisconnect = posToDisconnect.iter;

        if (nodeToDisconnect != head) {
            nodeToDisconnect->previous->next = nodeToDisconnect->next;
            nodeToDisconnect->next->previous = nodeToDisconnect->prev;
        }
    }

    bool eraseId(std::string id) {
        for(Node* currentNode(head->next); currentNode != head; currentNode = currentNode->next){
            if(currentNode->data.id == id) {
                currentNode->previous->next = currentNode->next;
                currentNode->next->previous = currentNode->previous;
                delete currentNode;
                --numberOfElements;
                return true;
            }
        }
        return false;
    }

    iterator erase(const iterator& posToDelete) {
        iterator result = posToDelete;
        ++result;

        Node* nodeToDelete = posToDelete.iter;

        if (nodeToDelete != head) {

            nodeToDelete->previous->next = nodeToDelete->next;
            nodeToDelete->next->previous = nodeToDelete->previous;

            delete nodeToDelete;
            --numberOfElements;
        }
        return result;
    }
    iterator erase(iterator& first, const iterator& last) {
        iterator result{ end() };
        if (first == begin() && last == end())
            clear();
        else {
            while (first != last)
                first = erase(first);
            result = last;
        }
        return result;
    }

    void trim_left(size_t n){
        size_t counter = 0;
        for(iterator it = begin(); it != end(); ++it){
            if(counter++ < n){
                erase(it);
            }
        }
    }

    List<T> copy_part(iterator& first, const iterator& last) {
        List<T> part;
        while (first != last) {
            part.push_back(first.iter->data);
            first++;
        }
        return part;
    }

    List<T> copy_part(int index_first, int index_last) {
        return copy_part(begin() + index_first, begin() + index_last);
    }

    List<T> grab_part(iterator& first, const iterator& last) {
        List<T> part;
        while (first != last) {
            part.push_back(first.iter->data);
            first = erase(first);
        }
        return part;
    }

    void pop_front() { erase(begin()); };
    void push_front(const T& d) { insert(begin(), d); }

    void pop_back() { erase(--end()); };
    void push_back(const T& d) { insert(end(), d); }

    void emplace_back(List<T>& other, iterator n) {
        n.iter->next->previous = n.iter->previous;
        n.iter->previous->next = n.iter->next;
        
        n.iter->previous = head->previous;
        n.iter->next = head;

        head->previous->next = n.iter;
        head->previous = n.iter;

        other.numberOfElements--;
        numberOfElements++;
    }
    
    void resize(size_t count, const T& value) {
        if (numberOfElements < count)
            for (size_t i{ numberOfElements }; i < count; ++i)
                insert(end(), value);
        else
            while (count--)
                pop_back();
    }
    void resize(size_t count) {
        if (numberOfElements < count)
            for (size_t i{ numberOfElements }; i < count; ++i)
                insert(end());
        else
            while (count--)
                pop_back();
    }

    void append(const List<T>& li) {
        for (List<T>::iterator it = li.begin(); it != li.end(); ++it) {
            push_back(it.iter->data);
        }
    }

    void append(List&& other){
        if(other.empty()) return;

        numberOfElements += other.numberOfElements;
        head->previous->next = other.head->next;
        other.head->next->previous = head->previous;

        head->previous = other.head->previous;
        other.head->previous->next = head;

        other.head->next = nullptr;
        other.head->previous = nullptr;
        delete other.head;

        other.init();
    }

    void swap(List& other) { std::swap(head, other.head); std::swap(numberOfElements, other.numberOfElements); }

    // Operations --------------------
    void reverse() {
        const Node* oldHead = head;

        for (Node* nptr = head; ; nptr = nptr->previous) {
            std::swap(nptr->next, nptr->previous);
            if (nptr->previous == oldHead) // Previous was the original next
                break;
        }
    }

    std::vector<T> tovec(){
        std::vector<T> v;
        for(Node* curr = head->next; curr != head; curr=curr->next){
            v.push_back(*curr);
        }
    }

    template<typename Lambda>
    void foreach(Lambda func) { // or Lambda&&, which is usually better
        Node* curr = head;
        while (curr != nullptr) {
            func(curr);
            curr = curr->next;
        }
    }

    // Non standard inefficient functions --------------------------
    T& operator[](const size_t index) const { return begin()[index]; }
    iterator at(const size_t index) const { return begin() + index; }

    struct const_iterator {

        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = T;
        using pointer = T*;
        using reference = T&;

        const Node* iter{};
        const Node* head{};

    };

    // ------------------------------------------------------------------------
    // Define iterator capability ---------------------------------------------
    struct iterator {

        // Definitions ----------------
        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = T;
        using pointer = T*;
        using reference = T&;

        // Data -----------------------
        Node* iter{};
        Node* head{};

        // Constructor ----------------
        iterator(Node* const node, Node* const h) : iter(node), head(h) {};
        iterator() {};

        // Dereferencing --------------
        reference operator*() const { return iter->data; }
        reference operator->() const { return &**this; }

        // Arithmetic operations ------
        iterator operator++() { iter = iter->next; return *this; }
        iterator operator--() { iter = iter->previous; return *this; }
        iterator operator++(int) { iterator tmp = *this; ++* this; return tmp; }
        iterator operator--(int) { iterator tmp = *this; --* this; return tmp; }

        iterator next() { iterator temp { *this }; return ++temp; }
        iterator prev() { iterator temp { *this }; return --temp; }

        iterator operator +(const difference_type& n) const {
            iterator temp{ *this };  difference_type k{ n }; if (k > 0) while (k--)++temp; else while (k++)--temp; return temp;
        }
        iterator operator +=(const difference_type& n) {
            difference_type k{ n }; if (k > 0) while (k--)++* this; else while (k++)--* this; return *this;
        };
        iterator operator -(const difference_type& n) const {
            iterator temp{ *this };  difference_type k{ n }; if (k > 0) while (k--)--temp; else while (k++)++temp; return temp;
        }
        iterator operator -=(const difference_type& n) {
            difference_type k{ n }; if (k > 0) while (k--)--* this; else while (k++)++* this; return *this;
        };
        // Comparison ----------------- (typical space ship implementation)
        bool operator ==(const iterator& other) const { return iter == other.iter; };
        bool operator !=(const iterator& other) const { return iter != other.iter; };
        bool operator < (const iterator& other) const { return other.iter - iter < 0; };
        bool operator <= (const iterator& other) const { return other.iter - iter <= 0; };
        bool operator > (const iterator& other) const { return other.iter - iter > 0; };
        bool operator >= (const iterator& other) const { return other.iter - iter >= 0; };

        // Special non standard functions -----------------
        difference_type operator-(const iterator& other) const {
            difference_type result{};
            Node* nptr = head;

            int indexThis{ -1 }, indexOther{ -1 }, index{};

            do {
                nptr = nptr->next;
                if (nptr == iter)
                    indexThis = index;
                if (nptr == other.iter)
                    indexOther = index;
                ++index;
            } while (nptr != head);

            if (indexThis >= 0 and indexOther >= 0)
                result = indexThis - indexOther;
            return result;
        }
        reference operator[] (const size_t index) {
            Node* nptr = head->next;
            for (size_t i{}; i < index and nptr != head; ++i, nptr = nptr->next)
                ;
            return nptr->data;
        }
    };


};

//class CustomListTA : public CustomList<TA> {
//public:
//    CustomListTA() {}
//    CustomListTA(const std::vector<TA>& v) { init(); insert(begin(), v.begin(), v.end()); }
//    CustomListTA(const std::vector<TA*>& v) { init(); for (auto& ta : v) push_back(*ta); }
//    CustomListTA(const CustomList& other) { init(); insert(begin(), other.begin(), other.end()); }
//    CustomListTA(CustomList&& other) { init(); head = other.head; other.clear(); }
//    CustomListTA(const iterator& first, const iterator& last) { init(); insert(begin(), first, last); }
//
//    CustomListTA& operator =(const CustomList& other) { clear(); insert(begin(), other.begin(), other.end()); return *this; }
//    CustomListTA& operator =(CustomList&& other) { clear(); head = other.head; other.clear(); return *this; }
//
//    int collectProfit() const { int sum{}; for (auto n : *this) { sum += n.profit; } return sum; }
//
//} ;