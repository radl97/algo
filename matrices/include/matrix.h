#ifndef MATRIX_H
#define MATRIX_H

#include<iostream>
#include<initializer_list>
#include<array>
#include<cassert>
template<int N, int M, typename T=int>
class Matrix;

template<int SIZE, typename TYPE=int>
inline static Matrix<SIZE,SIZE,TYPE> Indentity();


template<int N, int M, typename T>
class Matrix
{
    private:

        std::array<T,N*M> matrix;
    public:


        Matrix(){}

        /**
        *Constructor for initialize with an initializer_list.
        *If the list contains 1 item, the whole matrix will be initialized with the value.
        *else the list should contain N*M item.
        */
        Matrix(const std::initializer_list<T> list)
        {
            if(list.size()==1)
            {
                for(int i=0;i<N*M;i++)
                {
                    matrix[i]=*list.begin();
                }
                return;
            }
            assert(list.size()==N*M);
            int i=0;
            for(T item:list)
            {
                matrix[i]=item;
                i++;
            }
        }

        ///GET ITEM FROM MATRIX
        T& operator() (const int i,const int j)
        {
            return matrix[i*M+j];
        }

        T operator() (const int i,const int j) const
        {
            return matrix[i*M+j];
        }

        T& operator() (const int i)
        {
            return matrix[i];
        }

        T operator() (const int i) const
        {
            return matrix[i];
        }

        ///GET DIMENSION OF MATRIX
        int getWitdth() const
        {
            return M;
        }

        int getHeigth() const
        {
            return N;
        }

        ///MATRIX OPERATIONS

        Matrix& operator+=(const Matrix<N,M,T>& m2)
        {
            for(int i=0;i<N*M;i++)
            {
                (*this)(i)+=m2(i);
            }
            return *this;
           // Matrix* m=new Matrix<N,M();
        }

        Matrix operator+(const Matrix<N,M,T>& m2) const
        {
            Matrix tmp=(*this);
            return tmp+=m2;
        }


        template<int K>
        Matrix operator *(const Matrix<M,K,T>& m2) const
        {
            Matrix<N,K,T> c={0};

            for (int i=0; i<N; i++)
                for (int k=0; k<M; k++)
                    for (int j=0; j<M; j++)
                        c(i,j) += this->operator()(i,k)*m2(k,j);
            return c;
        }

        /**
        *Concants two matrices into the return value
        */
        template<int K>
        Matrix<N,M+K,T> concat(const Matrix<N,K,T>& m2) const
        {
            Matrix<N,M+K,T> res;
            for(int i=0;i<N;i++)
            {
                for(int j=0;j<M;j++)
                    res(i,j)=(*this)(i,j);

                for(int j=0;j<K;j++)
                    res(i,M+j)=m2(i,j);
            }
            return res;
        }

        /**
        *Return a new matrix with the last K column of this matrix
        */
        template<int K>
        Matrix<N,K,T> lastNColumn()
        {
            Matrix<N,K,T> m;
            for(int i=0;i<N;i++)
            {
                for(int j=0;j<M;j++)
                {
                    m(i,j)=(*this)(i,M-K+j);
                }
            }
            return m;
        }

        ///GAUSSIAN-ELIMINATION
        /**
        *Gaussian elimination using augmented form.
        *@warning Highly advised to use on double matrix
        *@return true if elimination was successful
        */
        bool EliminateAugmented(int augmented=1)
        {

            if(N<M-augmented)
                return false;
            //FIST PHASE
            for(int i=0;i<N && i<M-augmented;i++) //i: current row, stops when reaches last column as well
            {
                int r=firstNotNull(i);
                if(r==-1) //better with small DELTA_epsilon?
                    return false;

                //std::cout<<"before swap"<<std::endl;
                //print();
                T first=(*this)(r,i);
                for(int j=i;j<M;j++) //swaps the current row with the max row
                {
                    T tmp=(*this)(r,j)/first; //devide row as well
                    (*this)(r,j)=(*this)(i,j);
                    (*this)(i,j)=tmp;
                }
                //std::cout<<"after swap"<<std::endl;
                //print();

                //subtract current row from rows below;
                for(int k=i+1;k<N;k++) //k:row from which i. row is subtracted
                {
                    T c=-(*this)(k,i)/(*this)(i,i);
                    for(int j=i;j<M;j++)
                    {
                        (*this)(k,j)+=c*(*this)(i,j);
                    }
                }
            }

            int hnn=firstNotNull(M-augmented);
            if(hnn!=-1) //remained row with nonzero echelon
                return false;

            //SECOND PHASE
            for(int i=M-augmented-1;i>=0;i--) //starts from last column
            {
                for(int j=i-1;j>=0;j--)
                {
                    for(int k=M-augmented;k<M;k++)
                    {
                        (*this)(j,k)-=(*this)(j,i)*(*this)(i,k);
                    }
                    (*this)(j,i)=0;

                }
            }

            return true;
        }

        /**
        *Eliminate a square matrix.
        *@return determinant of the matrix
        */
        int EliminateSquare()
        {
            assert(N==M);
            int determinant=1;
            for(int i=0;i<N;i++)
            {
                int r=firstNotNull(i);
                if(r==-1)
                    return 0;
                if(r!=i)
                    determinant*=-1;

                T first=(*this)(r,i);
                determinant*=first;
                for(int j=i;j<N;j++) //swaps the current row with the max row
                {
                    T tmp=(*this)(r,j)/first; //devide row as well
                    (*this)(r,j)=(*this)(i,j);
                    (*this)(i,j)=tmp;
                }


                 for(int k=i+1;k<N;k++) //k:row from which i. row is subtracted
                {
                    T c=-(*this)(k,i)/(*this)(i,i);
                    for(int j=i;j<N;j++)
                    {
                        (*this)(k,j)+=c*(*this)(i,j);
                    }
                }

            }
            return determinant;

        }

        /**
        *Returns the row with the maximal value at the i. column.
        *Starts from the i. row
        */
        int firstNotNull(int i) const
        {
            int j;
            for(j=i;j<N && (*this)(j,i)==0;j++);
            return (j==N)?-1:j;
        }

        /**
        *Invert a matrix
        */
        void Invert()
        {
            assert(N==M);
            assert(det((*this))!=0);

            Matrix<N,2*N,T> m2=this->concat(Indentity<N,double>());
            m2.EliminateAugmented(N);
            (*this)=m2.template lastNColumn<N>();

        }



        ///DEBUG

        /**
        *Prints the matrix to the standard input with one extra line for separation
        */
        void print() const
        {
            for(int i=0;i<N;i++)
            {
                for(int j=0;j<M;j++)
                {
                    std::cout<<(*this)(i,j)<<" ";
                }
                std::cout<<std::endl;
            }
            std::cout<<std::endl;
        }

    protected:


};


/**
*Returns an indenity matrix with the given size, and given numeric type
*/
template <int SIZE, typename TYPE=int>
inline static Matrix<SIZE,SIZE,TYPE> Indentity()
{
    Matrix<SIZE,SIZE,TYPE> m{0};
    for(int i=0;i<SIZE;i++)
    {
        m(i,i)=1;
    }
    return m;
}

/**
*Returns a nullmatrix
*/
template<int N, int M, typename TYPE=int>
inline static Matrix<N,M,TYPE> NullMatrix()
{
    Matrix<N,M,TYPE> m={0};
    return m;
}

/**
*Calculates the determinant of a matrix
*@return the determinant of the matrix
*Copies the matrix, so it won't be altered.
*/
template<int N>
inline static int det(const Matrix<N,N,double>& m)
{
    Matrix<N,N,double> m2=m;
    return m2.EliminateSquare();
}

/**
*Calculates the inverse of a matrix, its content won't be altered
*@return inverse of the given matrix
*/
template<int N, typename T>
inline static Matrix<N,N> inverse(const Matrix<N,N,T>& m)
{
    assert(det(m)!=0);

    Matrix<N,2*N,T> m2=m->concat(Indentity<N,double>());
    m2.EliminateAugmented(N);
    return m2.template lastNColumn<N>();

}

#endif // MATRIX_H
