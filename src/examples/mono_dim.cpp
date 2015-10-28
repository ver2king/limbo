//#define SHOW_TIMER
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>

#include <limbo/limbo.hpp>
#include <limbo/inner_opt/cmaes.hpp>

using namespace limbo;

struct Params {
    struct gp_ucb : public defaults::gp_ucb {
    };

    struct cmaes : public defaults::cmaes {
    };

    struct ucb {
        BO_PARAM(float, alpha, 0.1);
    };

    struct kf_maternfivehalfs {
        BO_PARAM(float, sigma, 1);
        BO_PARAM(float, l, 0.2);
    };

    struct boptimizer {
        BO_PARAM(double, noise, 0.001);
        BO_PARAM(int, dump_period, 1);
    };

    struct init {
        BO_PARAM(int, nb_samples, 5);
    };

    struct maxiterations {
        BO_PARAM(int, n_iterations, 20);
    };
};

struct fit_eval {
    static constexpr size_t dim_in = 2;
    static constexpr size_t dim_out = 1;

    Eigen::VectorXd operator()(const Eigen::VectorXd& x) const
    {
        Eigen::VectorXd res(1);
        res(0) = 0;
        for (int i = 0; i < x.size(); i++)
            res(0) += 1 - (x[i] - 0.3) * (x[i] - 0.3) + sin(10 * x[i]) * 0.2;
        return res;
    }
};

int main()
{

    typedef kernel::MaternFiveHalfs<Params> Kernel_t;
    typedef mean::Data<Params> Mean_t;
    typedef model::GP<Params, Kernel_t, Mean_t> GP_t;
    typedef acqui::UCB<Params, GP_t> Acqui_t;

    bayes_opt::BOptimizer<Params, modelfun<GP_t>, acquifun<Acqui_t>> opt;
    opt.optimize(fit_eval());
    std::cout << opt.best_observation() << " res  "
              << opt.best_sample().transpose() << std::endl;
    return 0;
}
