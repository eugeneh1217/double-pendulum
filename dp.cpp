/**
 * double pendulum math
 * - https://web.mit.edu/jorloff/www/chaosTalk/double-pendulum/double-pendulum-en.html
 * 
 * double pendulum energy
 * - https://scienceworld.wolfram.com/physics/DoublePendulum.html
 * 
 * runge-kutta
 * - https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
*/

#define _USE_MATH_DEFINES
#define DEBUG

#include <SDL2/SDL.h>
#include <iostream>
#include <cmath>
#include <chrono>
#include <ctime>
#include <thread>
#include <functional>
#include <array>

// Simulation Constants
#define FRAME_RATE 20
#define GRAVITY 9.81
#define M1 10
#define M2 10
#define L1 1
#define L2 1
#define THETA1 M_PI/2 
#define THETA2 M_PI/2

// GUI Constants
const int SCREEN_WIDTH = 2160;
const int SCREEN_HEIGHT = 1440;
#define XSCALE 500
#define YSCALE 500
#define XOFFSET SCREEN_WIDTH / 2
#define YOFFSET SCREEN_HEIGHT * 1 / 5

// window
SDL_Window* gWindow = NULL;

// renderer
SDL_Renderer* gRenderer = NULL;

class Color
{
    public:
    unsigned r;
    unsigned g;
    unsigned b;
    unsigned a;

    Color(unsigned r, unsigned g, unsigned b, unsigned a)
    {
        this->r = r;
        this->g = g;
        this->b = b;
        this->a = a;
    }

    Color(unsigned r, unsigned g, unsigned b)
    {
        this->r = r;
        this->g = g;
        this->b = b;
        this->a = 0xFF;
    }
};

template<class T>
class Vector2
{
    public:
    T x;
    T y;

    Vector2(T x, T y)
    {
        this->x = x;
        this->y = y;
    }

    Vector2()
    {
        x = 0;
        y = 0;
    }
};

class Line
{
    public:
    Vector2<double> start;
    double angle;
    double length;
    double mass;
    double angular_speed;

    Line(Vector2<double> start, double angle, double length, double mass, double angular_speed)
    {
        this->start = start;
        this->angle = angle;
        this->length = length;
        this->mass = mass;
        this->angular_speed = angular_speed;
    }

    Line(Vector2<double> start, double angle, double length, double mass)
    {
        this->start = start;
        this->angle = angle;
        this->length = length;
        this->mass = mass;
        this->angular_speed = 0;
    }

    Line()
    {
        this->start = Vector2<double>();
        angle = 0;
        length = 0;
        mass = 0;
        angular_speed = 0;
    }

    void Rotate(double const &angle_delta)
    {
        this->angle = std::fmod((this->angle + angle_delta), (double) 2 * M_PI);
    }

    Vector2<double> GetEnd()
    {
        return Vector2<double>(
            start.x + length * std::cos(angle),
            start.y + length * std::sin(angle)
        );
    }
};

int compute_step(Line &l1, Line &l2, double deltatime)
{
    double t1 = l1.angle - 3.*M_PI/2;
    double t2 = 3.*M_PI/2 - l2.angle;
    double g = GRAVITY;
    double a1 = -g * (2 * l1.mass + l2.mass) * std::sin(t1);
    double b1 = l2.mass * g * std::sin(t1 - 2 * t2);
    double c1 = std::pow(l2.angular_speed, 2) * l2.length;
    double d1 = std::pow(l1.angular_speed, 2) * l1.length * std::cos(t1 - t2);
    double e1 = 2 * sin(t1 - t2) * l2.mass * (c1 + d1);
    double f = 2 * l1.mass + l2.mass - l2.mass * std::cos(2 * t1 - 2 * t2);
    double aa1 = (a1 - b1 - e1) / (l1.length * f);

    double a2 = 2 * std::sin(t1 - t2);
    double b2 = std::pow(l1.angular_speed, 2) * l1.length * (l1.mass + l2.mass);
    double c2 = g * (l1.mass + l2.mass) * std::cos(t1);
    double d2 = std::pow(l2.angular_speed, 2) * l2.length * l2.mass * std::cos(t1 - t2);
    double aa2 = a2 * (b2 + c2 + d2) / (l2.length * f);

    l1.angle += l1.angular_speed * deltatime;
    l2.angle += l2.angular_speed * deltatime;
    l1.angular_speed += aa1 * deltatime;
    l2.angular_speed += -aa2 * deltatime;

    l2.start = l1.GetEnd();

    return 0;
}

class double_pendulum
{
    public:
    double t;

    double_pendulum(Line l1, Line l2)
    {
        _l1 = l1;
        _l2 = l2;
        state[0] = l1.angle;
        state[1] = l1.angular_speed;
        state[2] = l2.angle;
        state[3] = l2.angular_speed;
    }

    void step(double step_size)
    {
        _runge_kutta(t, state, step_size, state);
        _l1.angle = state[0];
        _l1.angular_speed = state[1];
        _l2.angle = state[2];
        _l2.angular_speed = state[3];
        _l1.angle = std::fmod(_l1.angle - M_PI/2, 2*M_PI);
        _l2.start = _l1.GetEnd();
        _l2.angle = std::fmod(_l2.angle - M_PI/2, 2*M_PI);
        t += step_size;
    }

    Line get_line1()
    {
        return _l1;
    }

    Line get_line2()
    {
        return _l2;
    }

    private:
    std::array<double, 4> state; // [theta1, omega1, theta2, omega2]
    Line _l1;
    Line _l2;
    std::array<double, 4> _next_state;
    void _runge_kutta(
        double t,
        std::array<double, 4> const &x,
        double h,
        std::array<double, 4> &x_next);
    void _step(double, std::array<double, 4> const &, std::array<double, 4> &);
};

void double_pendulum::_runge_kutta(
        double t,
        std::array<double, 4> const &x,
        double h,
        std::array<double, 4> &x_next)
{
    std::array<double, 4> k1, x_temp, k2, k3, k4;
    _step(t, x, k1);

    for (unsigned long i = 0; i < 4; ++ i)
    {
        x_temp[i] = x[i] + h * k1[i] / 2;
    }
    _step(t + (h/2), x_temp, k2);

    for (unsigned long i = 0; i < 4; ++ i)
    {
        x_temp[i] = x[i] + h * k2[i] / 2;
    }
    _step(t + (h/2), x_temp, k3);

    for (unsigned long i = 0; i < 4; ++ i)
    {
        x_temp[i] = x[i] + h * k3[i];
    }
    _step(t + (h/2), x_temp, k4);

    for (unsigned long i = 0; i < 4; ++ i)
    {
        x_next[i] = x[i] + (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}

void double_pendulum::_step(
    double t,
    std::array<double, 4> const &x,
    std::array<double, 4> &x_d)
{
    double t1, w1, t2, w2;
    t1 = x[0];
    w1 = x[1];
    t2 = x[2];
    w2 = x[3];

    x_d[0] = w1;
    x_d[2] = w2;

    double g = GRAVITY;
    double a1 = -g * (2 * _l1.mass + _l2.mass) * std::sin(t1);
    double b1 = -_l2.mass * g * std::sin(t1 - 2 * t2);
    double c1 = std::pow(w2, 2) * _l2.length;
    double d1 = std::pow(w1, 2) * _l1.length * std::cos(t1 - t2);
    double e1 = -2 * sin(t1 - t2) * _l2.mass * (c1 + d1);
    double f = 2 * _l1.mass + _l2.mass - _l2.mass * std::cos(2 * t1 - 2 * t2);
    x_d[1] = (a1 + b1 + e1) / (_l1.length * f);

    double a2 = 2 * std::sin(t1 - t2);
    double b2 = std::pow(w1, 2) * _l1.length * (_l1.mass + _l2.mass);
    double c2 = g * (_l1.mass + _l2.mass) * std::cos(t1);
    double d2 = std::pow(w2, 2) * _l2.length * _l2.mass * std::cos(t1 - t2);
    x_d[3] = a2 * (b2 + c2 + d2) / (_l2.length * f);
}

bool init()
{
    if (SDL_Init( SDL_INIT_VIDEO ) < 0)
    {
        std::cout << "SDL could not initialize: " << SDL_GetError() << std::endl;
        return false;
    }

    // create window
    gWindow = SDL_CreateWindow( "Double Pendulum", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN );
    if( gWindow == NULL )
    {
        std::cout << "Window could not be created! SDL_Error: " << SDL_GetError() << std::endl;
        return false;
    }

    gRenderer = SDL_CreateRenderer( gWindow, -1, SDL_RENDERER_ACCELERATED );
    if( gRenderer == NULL )
    {
        std::cout << "Renderer could not be created: " << SDL_GetError() << std::endl;
        return false;
    }

    // initialize renderer color white
    SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF);

    return true;
}

void close()
{
    // destroy window
    SDL_DestroyWindow( gWindow );
    gWindow = NULL;

    // quit SDL subsystems
    SDL_Quit();
}

Vector2<double> _physics_to_display_coordinates(Vector2<double> const &physics)
{
    return Vector2<double> ((XSCALE * physics.x) + XOFFSET, YOFFSET - (YSCALE * (physics.y)));
}

int draw_dp(SDL_Renderer *r, double_pendulum &dp, Color &c1, Color &c2)
{
    // adjust angles
    Line l1 = dp.get_line1();
    Line l2 = dp.get_line2();
    // draw upper pendulum
    Vector2<double> s1 = _physics_to_display_coordinates(l1.start);
    Vector2<double> e1 = _physics_to_display_coordinates(l1.GetEnd());
    SDL_SetRenderDrawColor(r, c1.r, c1.g, c1.b, c1.a);
    SDL_RenderDrawLine(r, (int) s1.x, (int) s1.y, (int) e1.x, (int) e1.y);
    // draw lower pendulum
    Vector2<double> s2 = _physics_to_display_coordinates(l2.start);
    Vector2<double> e2 = _physics_to_display_coordinates(l2.GetEnd());
    SDL_SetRenderDrawColor(r, c2.r, c2.g, c2.b, c2.a);
    SDL_RenderDrawLine(r, (int) s2.x, (int) s2.y, (int) e2.x, (int) e2.y);
}

int main()
{
    if (!init())
        return 1;
    // main loop flag
    bool quit = false;

    // event handler
    SDL_Event e;

    // physics objects
    Line l1(Vector2<double>(0, 0), THETA1, L1, M1);
    Color l1c(77, 96, 128);
    Line l2(l1.GetEnd(), THETA1, L2, M2);
    Color l2c(251, 149, 60);

    double_pendulum dp(l1, l2);

    // deltatime tracking
    auto last_time = std::chrono::high_resolution_clock::now();
    auto this_time = std::chrono::high_resolution_clock::now();
    float deltatime = 0;

    // frame count
    unsigned long frame_n = 0;

    // while application is running
    while(1)
    {
        // handle events on queue
        while( SDL_PollEvent( &e ) != 0 )
        {
            //user requests quit
            if( e.type == SDL_QUIT )
            {
                quit = true;
            }
        }

        if (quit)
            break;

        // get deltatime
        this_time = std::chrono::high_resolution_clock::now();
        deltatime = ((double) std::chrono::duration_cast<std::chrono::milliseconds>(this_time - last_time).count()) / 1000.;
        last_time = std::chrono::high_resolution_clock::now();

        #ifdef DEBUG
        std::cout << "frame: " << frame_n;
        std::cout << ", Total Energy: " <<
        0.5 * l1.mass * l1.length * l1.length * l1.angular_speed * l1.angular_speed // kinetic energy of mass 1 = 1/2*m*r^2*w^2
        + 0.5 * l2.mass * (
            l1.length * l1.length * l1.angular_speed * l1.angular_speed
            + l2.length * l2.length * l2.angular_speed * l2.angular_speed
            + 2 * l1.length * l2.length * l1.angular_speed * l2.angular_speed * std::cos(l1.angle - l2.angle)
        ) // kinetic energy of mass 2
        + l1.mass * GRAVITY * (l1.length + l2.length + std::sin(l1.angle)) // GPE of mass 1
        + l2.mass * GRAVITY * (l1.length + l2.length + std::sin(l1.angle) + std::sin(l2.angle)) // GPE of mass 2
        << " Joules";
        std::cout << ", deltatime: " << deltatime;
        std::cout << ", height 1: " << l1.mass * GRAVITY * (l1.length + l2.length + std::sin(l1.angle));
        std::cout << ", height 2: " << l2.mass * GRAVITY * (l1.length + l2.length + std::sin(l1.angle) + std::sin(l2.angle));
        std::cout << std::endl;
        #endif

        // clear screen
        SDL_SetRenderDrawColor( gRenderer, 33, 37, 41, 0xFF);
        SDL_RenderClear( gRenderer );

        // physics
        dp.step(deltatime);
        l1 = dp.get_line1();
        l2 = dp.get_line2();

        // render
        draw_dp(gRenderer, dp, l1c, l2c);

        // update screen
        SDL_RenderPresent( gRenderer );

        // framerate control
        std::this_thread::sleep_for(std::chrono::milliseconds((int) (1000./FRAME_RATE)));
        frame_n ++;
    }

    close();

    return 0;
}
