#pragma once

#include <type_traits>

namespace GCLF
{
/* SFINAE */
/* Substitution failure is not an error. */

// Some useful websites:
// 1. How to enable/disable member functions
// https://www.foonathan.net/2016/12/conditionally-removing-functions/

/// @brief  Check if a class has certain function.
/// @note Only when "func" is unique (not overloaded), use this method.
template <typename T>
class has_func
{
  // here "&C::func" is a pointer to function.
  template <typename C> static std::true_type test(decltype(&C::func));
  template <typename C> static std::false_type test(...);
public:
  static const bool value = decltype(test<T>(nullptr))::value;
  // If T has member function "func", test<T>(0) matches the first definition.
  // Otherwise, test<T>(0) matches the second definition.
  // The failure of the fisrt definition is not an error.
};

/// @brief  Check if a class has certain function.
/// @tparam Args Parameter types passed to the function.
/// @note Only when "func" is not unique (overloaded), use this method.
template <typename T, typename... Args>
class has_overloaded_func
{
  template <typename C, typename = decltype(std::declval<C>().func(std::declval<Args>()...))> static std::true_type test(int);
  template <typename C> static std::false_type test(...);
public:
  static constexpr bool value = decltype(test<T>(0))::value;
};

/// @brief  check if a class has certain function.
template <typename T, typename R>
class has_static_func
{
  // A non-static member function receives an instance of the class as first parameter.
  // So we need to explicitly declare the function type of "func".
  // "typename F" is the type of function, "F" is the instance of function.
  template<typename F, F> helper;
  // find the func. when use this, change the template parameters of helper.
  template <typename C> static constexpr R test(helper<void (*)(void), C::func>*);
  // return a default value, indicating failure.
  template <typename C> static constexpr R test(...);
public:
  static constexpr R value = test<T>(nullptr);
};

/// @brief  check if a class has certain function and get its return value.
template <typename T, typename R>
class has_static_func_get_return_val
{
  // A non-static member function receives an instance of the class as first parameter.
  // So we need to explicitly declare the function type of "func".
  // "typename F" is the type of function, "F" is the instance of function.
  template<typename F, F> helper;
  // find the func and get its return value. when use this, change the template parameters of helper.
  template <typename C> static constexpr R test(helper<R (*)(void), C::func>*) { return C::func(); }
  // return a default value, indicating failure.
  template <typename C> static constexpr R test(...) { return R(0); }
public:
  static constexpr R value = test<T>(nullptr);
};

/// @brief  check if a class has certain type.
template <typename T>
class has_type
{
  template <typename C> static constexpr std::true_type test(typename C::value_type*);
  template <typename C> static constexpr std::false_type test(...);
public:
  static constexpr bool value = decltype(test<T>(nullptr))::value;
};

/// @brief  check if a class has certain type and get the type.
template <typename T>
class has_type_then_get_it
{
  template <typename C> static constexpr auto test(typename C::value_type*) -> typename C::value_type;
  template <typename C> static constexpr void test(...);
public:
  typedef decltype(test<T>(nullptr)) value_type;
};

}// namespace GCLF