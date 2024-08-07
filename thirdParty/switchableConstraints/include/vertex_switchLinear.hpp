#pragma once

#include "g2o/core/base_vertex.h"
#include <math.h>

class VertexSwitchLinear : public g2o::BaseVertex<1, double>
{

public:
    VertexSwitchLinear();

    virtual void setToOriginImpl();

    virtual void oplusImpl(const number_t* update);

    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
    virtual bool setEstimateDataImpl(const number_t* et);


    double x() const { return _x; };


    //! The gradient at the current estimate is always 1;
    double gradient() const { return 1; } ;

private:
    
    double _x;

};
