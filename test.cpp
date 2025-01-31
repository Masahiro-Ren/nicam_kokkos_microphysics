#include <iostream>
#include <tuple>
#include <type_traits>

template <typename... Types>
struct type_array {
    using as_tuple = std::tuple<Types...>;

    template <std::size_t I>
    using get = std::tuple_element_t<I, as_tuple>;

    static constexpr std::size_t size = sizeof...(Types); 
};

class Tag1 {};
class Tag2 {};

using types = type_array<Tag1, Tag2>;

template <class TagName>
struct functor {
public:
    void operator()() { std::cout << "default \n"; }
private:
    TagName tag;
};

template <>
struct functor<Tag1> {
public:
    void operator()() { std::cout << "tag1 \n"; }
private:
    Tag1 tag;
};

template <>
struct functor<Tag2> {
public:
    void operator()() { std::cout << "tag2 \n"; }
private:
    Tag2 tag;
};

int main()
{
    functor<types::get<1>>()();
    return 0;
}