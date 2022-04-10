#pragma once

#include <iterator>
#include <initializer_list>
#include <algorithm>
#include <iostream>
#include <type_traits>
#include <vector>
#include <string>

#define DEFAULT_ID "d"

//struct Point {
//    double lat;
//    double lon;
//
//    Point() {}
//    Point(double lat, double lon) : lat(lat), lon(lon) {}
//};
//
//typedef struct TouristAttraction {
//    std::string id;
//    Point p;
//    TouristAttraction() : p(0, 0) {}
//    TouristAttraction(std::string id, Point p) : id(id), p(p) {}
//} TA;

// ------------------------------------------------------------------------------------------------
// This would be in a header file -----------------------------------------------------------------

// Type trait helper to identify iterators --------------------------------------------------------
template<typename T, typename = void>
struct is_iterator { static constexpr bool value = false; };
template<typename T>
struct is_iterator<T, typename std::enable_if<!std::is_same<typename std::iterator_traits<T>::value_type, void>::value>::type> {
    static constexpr bool value = true;
};

// The CustomList class ---------------------------------------------------------------------------------
template <typename T>
class CustomList {
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
    CustomList() { init(); }
    explicit CustomList(const size_t count) { init(); insert(begin(), count); }
    explicit CustomList(const size_t count, const T& value) { init(); insert(begin(), count, value); };
    explicit CustomList(const std::vector<T> v) { init(); insert(begin(), v.begin(), v.end()); }
    explicit CustomList(const std::vector<T*> v) { init(); for (auto& i : v) push_back(*i); }
    template <typename Iter>
    CustomList(const Iter& first, const Iter& last) { init(); insert(begin(), first, last); }
    CustomList(const CustomList& other) { init(), insert(begin(), other.begin(), other.end()); };

    CustomList(CustomList&& other) : head(other.head), numberOfElements(other.numberOfElements) { other.init(); }
    CustomList(const std::initializer_list<T>& il) { init(); insert(begin(), il.begin(), il.end()); }
    template <int N> CustomList(T(&other)[N]) { init(); insert(begin(), std::begin(other), std::end(other)); }
    template <int N> CustomList(const T(&other)[N]) { init(); insert(begin(), std::begin(other), std::end(other)); }


    // Assignment ---------------------
    CustomList& operator =(const CustomList& other) { clear(); insert(begin(), other.begin(), other.end()); return *this; }
    CustomList& operator =(CustomList&& other) { clear(); head = other.head; numberOfElements = other.numberOfElements; other.init(); return *this; }
    CustomList& operator =(const std::initializer_list<T>& il) { clear(); insert(begin(), il.begin(), il.end()); return *this; }
    template <int N> CustomList& operator =(const T(&other)[N]) { clear(); insert(begin(), std::begin(other), std::end(other)); return *this; }
    template <int N> CustomList& operator =(T(&other)[N]) { clear(); insert(begin(), std::begin(other), std::end(other)); return *this; }

    template <typename Iter> void assign(const Iter& first, const Iter& last) { clear(); insert(begin(), first, last); }
    template <int N> void assign(const T(&other)[N]) { clear(); insert(begin(), std::begin(other), std::end(other)); }
    template <int N> void assign(T(&other)[N]) { clear(); insert(begin(), std::begin(other), std::end(other)); }
    void assign(const size_t count, const T& value) { clear(); insert(begin(), count, value); }
    void assign(const std::initializer_list<T>& il) { clear(); insert(begin(), il.begin(), il.end()); }

    // Destructor ---------------------
    ~CustomList() { clear(); }

    // Element Access -----------------
    T& front() { return *begin(); }
    T& back() { return *(--end()); }

    // Iterators ----------------------
    iterator begin() const { return iterator(head->next, head); }
    iterator end() const { return iterator(head, head); }

    // Capacity -----------------------
    inline size_t size() const { return numberOfElements; }
    inline bool empty() { return size() == 0; }

    // Modifiers ----------------------
    void clear() {
        for (Node* nextNode{}, * currentNode(head->next); currentNode != head; currentNode = nextNode) {
            nextNode = currentNode->next;
            delete currentNode;
        }
        init();
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



    CustomList<T> copy_part(iterator& first, const iterator& last) {
        CustomList<T> part;
        while (first != last) {
            part.push_back(first.iter->data);
            first++;
        }
        return part;
    }

    CustomList<T> copy_part(int index_first, int index_last) {
        return copy_part(begin() + index_first, begin() + index_last);
    }

    CustomList<T> grab_part(iterator& first, const iterator& last) {
        CustomList<T> part;
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

    //this can be improved by implementing a vector's emplace_back function
    void append(const CustomList<T>& li) {
        for (CustomList<T>::iterator it = li.begin(); it != li.end(); ++it) {
            push_back(it.iter->data);
        }
    }

    void swap(CustomList& other) { std::swap(head, other.head); std::swap(numberOfElements, other.numberOfElements); }

    // Operations --------------------
    void reverse() {
        const Node* oldHead = head;

        for (Node* nptr = head; ; nptr = nptr->previous) {
            std::swap(nptr->next, nptr->previous);
            if (nptr->previous == oldHead) // Previous was the original next
                break;
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
        //reference operator->() const { return &**this; }

        // Arithmetic operations ------
        iterator operator++() { iter = iter->next; return *this; }
        iterator operator--() { iter = iter->previous; return *this; }
        iterator operator++(int) { iterator tmp = *this; ++* this; return tmp; }
        iterator operator--(int) { iterator tmp = *this; --* this; return tmp; }

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