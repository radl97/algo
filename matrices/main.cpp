#include <iostream>
#include "matrix.h"

using namespace std;

int main()
{
    Matrix<2,2,int> m={2,1,3,4};
    Matrix<2,2,int> m2={2,1,2,1};
    Matrix<2,2,int> m4={0};
    auto m5=Indentity<2,int>();
    auto m3=m+m2;
    auto m6=m*m2;
    m3.print();
    m4.print();
    m5.print();
    m6.print();

    Matrix<3,3,double> m7={1, 2, 4,
                            1, 2, 4,
                           2, -4, 5};
    m7.print();
    cout<<(m7.EliminateAugmented()?"success":"error")<<std::endl;
    m7.print();

    Matrix<3,3,double> m8={1, 2, 4,
                            1, 1, 4,
                           2, -4, 5};
    m8.print();
    cout<<det(m8)<<std::endl;
    m8.print();

    m8.Invert();
    m8.print();
    return 0;
}
