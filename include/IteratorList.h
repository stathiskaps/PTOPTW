#include <iterator>
#include <initializer_list>
#include <algorithm>
#include <iostream>
#include <type_traits>
#include <vector>

template<typename T, typename = void>
struct is_iterator { static constexpr bool value = false; };
template<typename T>
struct is_iterator<T, typename std::enable_if<!std::is_same<typename std::iterator_traits<T>::value_type, void>::value>::type> {
    static constexpr bool value = true;
};

//the list class
template <typename T>
class List {
    //subclass for a node
    struct Node{
        T data{};
        Node* next{};
        Node* previous{};
        Node() {}
        Node(Node* const n, Node* const p):  next(n), previous(p) {}
        Node(Node* const n, Node* const p, const T& d) : next(n), previous(p), data(d) {}
    };

    //private list data and funcctions
    Node* head{};
    size_t numberOfElements{};
    void init() { head = new Node(); head->next = head; head->prev = head; numberOfElements = 0; }

public:
    struct iterator; //Forward declaration

    //Constructor ------------
    List() { init(); }
    explicit List(const size_t count) { init(); insert(begin(), count); }
    explicit List(const size_t count, const T& value) { init(); insert(begin(), count, value); };
    template <typename Iter>
    List(const Iter& first, const Iter& last) { init(); insert(begin(), first, last); }
    List(const List& other) { init(), insert(begin(), other.begin(), other.end()); }

    List(List&& other) : head(other.head), numberOfElements(other.numberOfElements) { other.init(); }
    List(const std::initializer_list<T>& il) { init(); insert(begin(), il.begin(), il.end()); }
    template <int N> List(T(&other)[N]) { init(); insert(begin(), std::begin(other), std::end())}
};