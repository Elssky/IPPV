#ifndef LhCDSCVX_EDGEFN_H
#define LhCDSCVX_EDGEFN_H


class EdgeFN {
public:
    int from, to, index;
    double cap, flow;
    EdgeFN(int from, int to, double cap, double flow, double index);
};


#endif //LhCDSCVX_EDGEFN_H
