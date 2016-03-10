#include <iostream>
#include <polymul.h>

int main(void)
{
    // Two third-degree polynomials of two variables
    polymul::polynomial<int, 2, 3> p1, p2;

    // One sixth-degree polynomials of two variables
    polymul::polynomial<int, 2, 6> p3;

    // Set p1 = 7*x^2y and p2 = 5*y^3
    int exp1[2] = {2, 1};
    int exp2[2] = {0, 3};
    p1[p1.term_index(exp1)] = 7;
    p2[p2.term_index(exp2)] = 5;

    // Multiply p1 and p2, put the result in p3
    polymul::mul(p3, p1, p2);

    // Look what the resulting x^2y^4 term is
    int exp3[2] = {2, 4};
    std::cout << "The x^2y^4 term has coefficient " << p3[p3.term_index(exp3)] << std::endl;

    // Print all coefficients
    std::cout << "All coefficients of p3:" << std::endl;
    polymul::print(p3);

    return 0;
}
