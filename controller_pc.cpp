#ifndef CONTROLLER_CPP 
#define CONTROLLER_CPP 

#include <math.h> 
#include <vector>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

typedef double (*read_func_t)();
typedef void (*write_func_t)(double);

// controller base class 
struct Controller
{
    read_func_t  reader;
    write_func_t  writer;
    Controller(read_func_t reader, write_func_t writer) : reader(reader), writer(writer) {}
    virtual void update() = 0;
};

// first order IIR controller base class  
struct IIRBaseController {
    double le, lo; // last error and last output 
    IIRBaseController() : le{}, lo{} {}
    virtual double transfer(double const err) = 0;
};

// first order IIR controller with zero and pole 
struct IIRFirstOrderController : public IIRBaseController {
    double const z, p; 
    IIRFirstOrderController(double const z, double const p) : z(z), p(p) {} 
    double transfer(double const e) override {
        double ret = ((1 + M_PI * p) * lo - (1 + M_PI * z) * le + (1 - M_PI * z) * e) / (1 - M_PI * p); 
        lo = ret; 
        le = e; 
        return ret; 
    }
};

// first order IIR controller with only pole 
struct IIRSinglePoleController : public IIRBaseController {
    double const p; 
    IIRSinglePoleController(double const p) : p(p) {} 
    double transfer(double const e) override {
        double ret = ((1 + M_PI * p) * lo + M_PI * (le + e)) / (1 - M_PI * p);
        lo = ret; 
        le = e; 
        return ret; 
    }
};

// IIR controller in cascade form 
template<int len_zeroes, int len_poles> 
struct IIRCascadeController : public Controller {
    static_assert(len_zeroes < len_poles, "There should be more poles than zeroes. ");

    double const overall_gain, lower, upper;
    double last_out;
    std::vector<IIRBaseController*> controllers;

    IIRCascadeController(read_func_t reader, write_func_t writer, double const (&zeroes)[len_zeroes], double const (&poles)[len_poles], double const overall_gain, double const lower, double const upper) : Controller(reader, writer), overall_gain(overall_gain), lower(lower), upper(upper), last_out(0.)
    {
        int i = 0; 
        for (;i < len_zeroes; ++i)
            controllers.push_back(new IIRFirstOrderController(zeroes[i], poles[i]));
        for (;i < len_poles; ++i)
            controllers.push_back(new IIRSinglePoleController(poles[i]));
    }

    void update() override {
        double new_out = reader(); // start with an error 
        for(auto p : controllers) 
            new_out = p -> transfer(new_out);            
        new_out = max(min(overall_gain * new_out, upper), lower);
        if (abs(new_out - last_out) > 1.) // lazy update 
        {
            writer(int(new_out));
            last_out = (int)new_out;
        }
    }
};

// wrapper for constructors 
template<int len_zeroes, int len_poles> 
IIRCascadeController<len_zeroes, len_poles> make_iir_cascade_controller(read_func_t reader, write_func_t writer, double const (&zeroes)[len_zeroes], double const (&poles)[len_poles], double const overall_gain, double const lower = -32768., double const upper = 32767.) {
    return IIRCascadeController<len_zeroes, len_poles>(reader, writer, zeroes, poles, overall_gain, lower, upper);
}

// deprecated PID controller 
#if 0 
struct PIDController : public Controller
{
    double kp, ki, kd;
    double last_error, error, integral;

    PIDController() = delete;
    PIDController(read_func_t reader, write_func_t writer, double kp, double ki, double kd) : Controller(reader, writer), kp(kp), ki(ki), kd(kd), last_error(), error(), integral()
    {
        this->reader = reader;
        this->writer = writer;
    }

    void update() override
    {
        error           = reader();
        int32_t dac_num = -(kp * error + ki * integral + kd * (error - last_error));
        if (dac_num > 32000)
            dac_num = 32000; // ensure single side, and anti-windup
        else if (dac_num < 0)
            dac_num = 0; // ensure single side, and anti-windup
        else
            integral += error;
        writer((double)dac_num);
        last_error = error;
    }
};
#endif 

// deprecated IIRcontroller in its polynomial form 
#if 0 
template <int err_size, int output_size>
struct IIRController : public Controller
{
    double ec[err_size];
    double oc[output_size];
    double const lower, upper;
    std::deque<double> errs, outputs;

    IIRController() = delete;
    IIRController(read_func_t reader, write_func_t writer, double const (&err_coef)[err_size], double const (&output_coef)[output_size], double const lower, double const upper) : Controller(reader, writer), lower(lower), upper(upper), errs(err_size, 0.), outputs(output_size, 0.)
    {
        memcpy(ec, err_coef, sizeof(err_coef[0]) * err_size);
        memcpy(oc, output_coef, sizeof(output_coef[0]) * output_size);
    }

    void update() override
    {
        double error = reader();
        errs.pop_front();
        errs.push_back(error);
        double new_out = 0.;
        for (int i = 0; i < err_size; ++i)
            new_out += errs[i] * ec[i];
        for (int i = 0; i < output_size; ++i)
            new_out += outputs[i] * oc[i];
        new_out = max(min(new_out, upper), lower);
        writer(new_out); 
        outputs.push_back(new_out);
        outputs.pop_front();
    }

    void update_sim() 
    {
        double error = reader() + 4.1e-4 * outputs[output_size-1];
        errs.pop_front();
        errs.push_back(error);
        double new_out = 0.;
        for (int i = 0; i < err_size; ++i)
            new_out += errs[i] * ec[i];
        for (int i = 0; i < output_size; ++i)
            new_out += outputs[i] * oc[i];
        Serial.printf("%le ", error);
        writer(new_out); 
        outputs.push_back(new_out);
        outputs.pop_front();
    }
};

template <int err_size, int output_size> 
IIRController<err_size, output_size> make_iir_controller(read_func_t reader, write_func_t writer, double const (&err_coef)[err_size], double const (&output_coef)[output_size], double const lower = -32768., double const upper = 32767.) {
    return IIRController<err_size, output_size>(reader, writer, err_coef, output_coef, lower, upper);
}
#endif 
#endif  // CONTROLLER_CPP