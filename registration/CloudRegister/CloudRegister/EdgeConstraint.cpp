#include "BaseType.h"
#include "EdgeConstraint.h"
#include "g2o/core/factory.h"
#include "g2o/stuff/macros.h"
namespace g2o
{

// G2O_REGISTER_TYPE(EDGE_PT_TO_PLANE_DIST,EdgePtToPlaneDist);


//read
bool EdgePtToPlaneDist::read(std::istream& is)  
{
    double p;
    is >> p;
    setMeasurement(p);

    for (int i = 0; i < 1; ++i)
    {
        for (int j = i; j < 1; ++j) 
        {
            is >> information()(i, j);
            if (i != j)
                information()(j, i) = information()(i, j);
        }
    }
    return true;
}

//write
bool EdgePtToPlaneDist::write(std::ostream& os) const 
{
    double p = measurement();
    os << p;
    for (int i = 0; i < 1; ++i)
      for (int j = i; j < 1; ++j)
        os << " " << information()(i, j);
    return os.good();
}


} //namespace g2o

