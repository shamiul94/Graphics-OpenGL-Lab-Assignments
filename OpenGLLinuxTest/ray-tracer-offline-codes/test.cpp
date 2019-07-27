//
// Created by shamiul93 on 7/22/19.
//

#include <bits/stdc++.h>

using namespace std;

class base {
public:
    int a, b;

    base() { a = b = 0; }

    base(int a, int b) {
        this->a = a;
        this->b = b;
    }

    virtual void draw() {
        cout << "a : " << a << " b : " << b << endl;
    }
};

class derived1 : public base {
public:
    int a, b, c;

    derived1(int a, int b, int c) {
        this->a = a;
        this->b = b;
        this->c = c;
    }

    void draw() {
        cout << "a : " << a << " b : " << b << " c : " << c << endl;
    }
};

template <typename T1, typename T2>
void myfunc(T1 x, T2 y, int d)
{
    x.draw();
    y.draw();
    cout << d << endl ;
}

void f(int &a){
    a = 9;
}

int main() {
    int x = 5;
    cout << x << endl ;
    f(x);
    cout << x << endl ;
}