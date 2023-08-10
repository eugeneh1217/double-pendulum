#define _USE_MATH_DEFINES
#define DEBUG

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <iostream>
#include <cmath>
#include <chrono>
#include <thread>
#include <functional>
#include <array>
#include <vector>
#include <string>

#define OK 0
#define ERR 1

// Simulation Constants
#define FRAME_RATE 20
#define GRAVITY 9.81
#define DP_N 1

// first pendulum initial conditions
#define M1 0.01
#define M2 0.01
#define L1 1
#define L2 1
#define THETA1 3*M_PI/2 + M_PI/2
#define THETA2 3*M_PI/2 + M_PI/2

// pendulum initial condition changes
#define DM1 0
#define DM2 0
#define DL1 0
#define DL2 0
#define DTHETA1 0.01
#define DTHETA2 0.01
#define DCR 3
#define DCG 1
#define DCB 2

// GUI Constants
const int SCREEN_WIDTH = 2160;
const int SCREEN_HEIGHT = 2160;
#define XSCALE 400
#define YSCALE 400
#define XOFFSET SCREEN_WIDTH / 2
#define YOFFSET SCREEN_HEIGHT * 2 / 5
#define FONT_PATH "media/Ubuntu-Th.ttf"
#define FONT_SIZE 80
#define ENERGY_DIGITS 7
// #define DISPLAY_ENERGY

class color
{
    public:
    unsigned r;
    unsigned g;
    unsigned b;
    unsigned a;

    color()
    {
        r = 0;
        g = 0;
        b = 0;
        a = 0;
    }

    color(unsigned r, unsigned g, unsigned b, unsigned a)
    {
        this->r = r % 256;
        this->g = g % 256;
        this->b = b % 256;
        this->a = a % 256;
    }

    color(unsigned r, unsigned g, unsigned b)
    {
        this->r = r % 256;
        this->g = g % 256;
        this->b = b % 256;
        this->a = 256;
    }
};

template<class T>
class vector2
{
    public:
    T x;
    T y;

    vector2(T x, T y)
    {
        this->x = x;
        this->y = y;
    }

    vector2()
    {
        x = 0;
        y = 0;
    }
};

class Line
{
    public:
    vector2<double> start;
    double angle;
    double length;
    double mass;
    double angular_speed;

    Line(vector2<double> start, double angle, double length, double mass, double angular_speed)
    {
        this->start = start;
        this->angle = angle;
        this->length = length;
        this->mass = mass;
        this->angular_speed = angular_speed;
    }

    Line(vector2<double> start, double angle, double length, double mass)
    {
        this->start = start;
        this->angle = angle;
        this->length = length;
        this->mass = mass;
        this->angular_speed = 0;
    }

    Line()
    {
        this->start = vector2<double>();
        angle = 0;
        length = 0;
        mass = 0;
        angular_speed = 0;
    }

    void Rotate(double const &angle_delta)
    {
        this->angle = std::fmod((this->angle + angle_delta), (double) 2 * M_PI);
    }

    vector2<double> GetEnd()
    {
        return vector2<double>(
            start.x + length * std::cos(angle),
            start.y + length * std::sin(angle)
        );
    }
};

class double_pendulum
{
    public:
    double t;

    double_pendulum()
    {
        _l1 = Line();
        _l2 = Line();
        state[0] = 0;
        state[1] = 0;
        state[2] = 0;
        state[3] = 0;
    }

    double_pendulum(Line l1, Line l2) // l1 and l2 in terms of UI angles
    {
        _l1 = l1;
        _l2 = l2;
        state[0] = l1.angle + M_PI/2; // convert to physics angles
        state[1] = l1.angular_speed;
        state[2] = l2.angle + M_PI/2; // convert to physics angles
        state[3] = l2.angular_speed;
    }

    void step(double step_size)
    {
        _runge_kutta(t, state, step_size, state);
        _l1.angle = state[0];
        _l1.angular_speed = state[1];
        _l2.angle = state[2];
        _l2.angular_speed = state[3];
        _l1.angle = std::fmod(_l1.angle - M_PI/2, 2*M_PI); // convert back to UI angles
        _l2.start = _l1.GetEnd();
        _l2.angle = std::fmod(_l2.angle - M_PI/2, 2*M_PI); // convert back to UI angles
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

    double get_ke()
    {
        return 0.5 * _l1.mass * _l1.length * _l1.length * _l1.angular_speed * _l1.angular_speed
            + 0.5 * _l2.mass * (
                _l1.length * _l1.length * _l1.angular_speed * _l1.angular_speed
                + _l2.length * _l2.length * _l2.angular_speed * _l2.angular_speed
                + 2 * _l1.length * _l2.length * _l1.angular_speed * _l2.angular_speed * std::cos(_l1.angle - _l2.angle)
            );
    }

    double get_ge()
    {
        return _l1.mass * GRAVITY * (_l1.length + _l2.length + std::sin(_l1.angle))
            + _l2.mass * GRAVITY * (_l1.length + _l2.length + std::sin(_l1.angle) + std::sin(_l2.angle));
    }

    std::array<double, 4> state; // [theta1, omega1, theta2, omega2]
    private:
    Line _l1;
    Line _l2;
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
    __attribute__((unused)) double t,
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

vector2<double> _physics_to_display_coordinates(vector2<double> const &physics)
{
    return vector2<double> ((XSCALE * physics.x) + XOFFSET, YOFFSET - (YSCALE * (physics.y)));
}

int draw_dp(SDL_Renderer *r, double_pendulum &dp, color &c1, color &c2)
{
    // adjust angles
    Line l1 = dp.get_line1();
    Line l2 = dp.get_line2();
    // draw upper pendulum
    vector2<double> s1 = _physics_to_display_coordinates(l1.start);
    vector2<double> e1 = _physics_to_display_coordinates(l1.GetEnd());
    SDL_SetRenderDrawColor(r, c1.r, c1.g, c1.b, c1.a);
    SDL_RenderDrawLine(r, (int) s1.x, (int) s1.y, (int) e1.x, (int) e1.y);
    // draw lower pendulum
    vector2<double> s2 = _physics_to_display_coordinates(l2.start);
    vector2<double> e2 = _physics_to_display_coordinates(l2.GetEnd());
    SDL_SetRenderDrawColor(r, c2.r, c2.g, c2.b, c2.a);
    SDL_RenderDrawLine(r, (int) s2.x, (int) s2.y, (int) e2.x, (int) e2.y);
    return 0;
}

class object
{
    public:
    double_pendulum physics;
    color c1;
    color c2;

    object(double_pendulum physics, color c1, color c2)
    {
        this->physics = physics;
        this->c1 = c1;
        this->c2 = c2;
    }

    object(double_pendulum physics, color c)
    {
        this->physics = physics;
        c1 = c;
        c2 = c;
    }

    object(double_pendulum physics)
    {
        this->physics = physics;
        c1 = color(77, 96, 128);
        c2 = c1;
    }
};

void print_sdl_error(std::string msg)
{
    std::cout << msg << ": " << SDL_GetError() << std::endl;
}

std::string double_to_string(double d, unsigned digits)
{
    // double r = std::round(d * std::pow(10, precision)) / std::pow(10, precision);
    // std::string s = std::to_string(r);
    // std::string clipped = s.substr(0, s.find(".") + precision + 1);
    // return clipped;
    return std::to_string(d).substr(0, digits + 1);
}

class sdl_simulation
{
    public:
    std::vector<object> dps;

    sdl_simulation()
    {
        dps = {};
        gWindow = NULL;
        gRenderer = NULL;
    }

    sdl_simulation(std::vector<object> dps)
    {
        this->dps = dps;
        gWindow = NULL;
        gRenderer = NULL;
    }

    int run()
    {
        if (_init() != OK)
            return 1;

        // main loop flag
        bool quit = false;

        // event handler
        SDL_Event e;

        // deltatime tracking
        auto last_time = std::chrono::high_resolution_clock::now();
        auto this_time = std::chrono::high_resolution_clock::now();
        float deltatime = 0;

        // frame count
        unsigned long frame_n = 0;

        #ifdef DISPLAY_ENERGY
        // text rendering
        SDL_Surface *text_surface;
        SDL_Texture *text_texture;
        std::string energy_str = "";
        SDL_Color tc = {255, 255, 255, 0};

        // energy
        double ge = 0;
        double ke = 0;
        #endif

        // while application is running
        while(1)
        {
            // handle events on queue
            while( SDL_PollEvent( &e ) != 0 )
            {
                //user requests quit
                if( e.type == SDL_QUIT )
                {
                    return 0;
                }
            }

            if (quit)
                break;
            
            // if (frame_n == 0)
            // {
            //     std::this_thread::sleep_for(std::chrono::seconds(15));
            //     last_time = std::chrono::high_resolution_clock::now();
            // }

            // get deltatime
            this_time = std::chrono::high_resolution_clock::now();
            deltatime = ((double) std::chrono::duration_cast<std::chrono::milliseconds>(this_time - last_time).count()) / 1000.;
            last_time = std::chrono::high_resolution_clock::now();

            #ifdef DEBUG
            std::cout << "frame: " << frame_n;
            std::cout << ", framerate: " << 1./deltatime;
            std::cout << std::endl;
            #endif

            // clear screen
            SDL_SetRenderDrawColor( gRenderer, 33, 37, 41, 0xFF);
            SDL_RenderClear( gRenderer );

            // update objects
            #ifdef DISPLAY_ENERGY
            ge = 0;
            ke = 0;
            #endif
            for (auto &dp : dps)
            {
                dp.physics.step(deltatime);
                if (_draw_dp(gRenderer, dp) != OK)
                {
                    quit = true;
                    break;
                }
                #ifdef DISPLAY_ENERGY
                ge += dp.physics.get_ge();
                ke += dp.physics.get_ke();
                #endif
            }

            #ifdef DISPLAY_ENERGY
            energy_str =
                "GPE: " + double_to_string(ge, ENERGY_DIGITS) + "J\n"
                + "KE: " + double_to_string(ke, ENERGY_DIGITS) + "J\n"
                + "Total Energy: " + double_to_string(ge + ke, ENERGY_DIGITS) + "J";
            text_surface = TTF_RenderText_Solid_Wrapped(font, energy_str.c_str(), tc, 2048);
            text_texture = SDL_CreateTextureFromSurface(gRenderer, text_surface);
            int tw = 0;
            int th = 0;
            SDL_QueryTexture(text_texture, NULL, NULL, &tw, &th);
            SDL_Rect text_rect = {0, 0, tw, th};
            SDL_RenderCopy(gRenderer, text_texture, NULL, &text_rect);
            #endif

            // update screen
            SDL_RenderPresent( gRenderer );

            #ifdef DISPLAY_ENERGY
            SDL_FreeSurface(text_surface);
            SDL_DestroyTexture(text_texture);
            #endif

            // frame control
            frame_n ++;
            std::this_thread::sleep_for(std::chrono::milliseconds((int) (1000./FRAME_RATE)));
        }

        _close();

        return OK;
    }

    private:
    SDL_Window* gWindow;
    SDL_Renderer* gRenderer;
    TTF_Font *font;

    int _init()
    {
        if (SDL_Init( SDL_INIT_VIDEO ) < 0)
        {
            print_sdl_error("SDL could not initialize");
            return ERR;
        }

        if (TTF_Init() < 0)
        {
            print_sdl_error("SDL tff could not initialize");
            return ERR;
        }

        // load font
        font = TTF_OpenFont(FONT_PATH, FONT_SIZE);

        // create window
        gWindow = SDL_CreateWindow( "Double Pendulum", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN );
        if( gWindow == NULL )
        {
            print_sdl_error("Window could not be created");
            return ERR;
        }

        gRenderer = SDL_CreateRenderer( gWindow, -1, SDL_RENDERER_ACCELERATED );
        if( gRenderer == NULL )
        {
            print_sdl_error("Renderer could not be created");
            return ERR;
        }

        // initialize renderer color white
        SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF);

        return OK;
    }

    int _draw_dp(SDL_Renderer *r, object obj)
    {
        Line l1 = obj.physics.get_line1();
        Line l2 = obj.physics.get_line2();
        int ret = 0;
        // draw upper pendulum
        vector2<double> s1 = _physics_to_display_coordinates(l1.start);
        vector2<double> e1 = _physics_to_display_coordinates(l1.GetEnd());
        ret = SDL_SetRenderDrawColor(r, obj.c1.r, obj.c1.g, obj.c1.b, obj.c1.a);
        if (ret != OK)
        {
            print_sdl_error("Could not switch color");
            return ERR;
        }
        ret = SDL_RenderDrawLine(r, (int) s1.x, (int) s1.y, (int) e1.x, (int) e1.y);
        if (ret != OK)
        {
            print_sdl_error("Could not draw line");
            return ERR;
        }
        // draw lower pendulum
        vector2<double> s2 = _physics_to_display_coordinates(l2.start);
        vector2<double> e2 = _physics_to_display_coordinates(l2.GetEnd());
        ret = SDL_SetRenderDrawColor(r, obj.c2.r, obj.c2.g, obj.c2.b, obj.c2.a);
        if (ret != OK)
        {
            print_sdl_error("Could not switch color");
            return ERR;
        }
        ret = SDL_RenderDrawLine(r, (int) s2.x, (int) s2.y, (int) e2.x, (int) e2.y);
        if (ret != OK)
        {
            print_sdl_error("Could not draw line");
            return ERR;
        }
        return OK;
    }

    void _close()
    {
        // destroy renderer
        SDL_DestroyRenderer(gRenderer);
        gRenderer = NULL;
        // destroy window
        SDL_DestroyWindow( gWindow );
        gWindow = NULL;

        // quit sdl_ttf
        TTF_CloseFont(font);
        TTF_Quit();

        // quit SDL subsystems
        SDL_Quit();
    }
};

int main()
{
    std::vector<object> objs{};
    Line l1 = Line(vector2<double> (0, 0), THETA1, L1, M1, 0);
    Line l2 = Line(l1.GetEnd(), THETA2, L2, M2, 0);
    color c1(77, 96, 128);
    // color c2(251, 149, 60);
    for (int i = 0; i < DP_N; ++ i)
    {
        l1.mass += DM1;
        l1.length += DL1;
        l1.angle += DTHETA1;
        l2.mass += DM2;
        l2.length += DL2;
        l2.angle += DTHETA2;
        c1.r += DCR;
        c1.g += DCG;
        c1.b += DCB;
        objs.push_back(object (double_pendulum(l1, l2), c1, c1));
    }
    sdl_simulation sim(objs);
    sim.run();

    return 0;
}
