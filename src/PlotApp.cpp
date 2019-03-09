#include <Magnum/Platform/Sdl2Application.h>

#include "FunctionPlot.h"

using namespace Magnum;
using namespace Magnum::Math::Literals;

class MyApp: public Platform::Application {
    public:
        explicit MyApp(const Arguments& arguments);

    private:
        void drawEvent() override;
        void mousePressEvent(MouseEvent& event) override;
        void mouseReleaseEvent(MouseEvent& event) override;
        void mouseMoveEvent(MouseMoveEvent& event) override;

        FunctionPlot _plot;
};

MyApp::MyApp(const Arguments& arguments):
    Platform::Application{arguments, Configuration{}.setTitle("Function Plot App"), GLConfiguration{}.setSampleCount(16)}, _plot{windowSize()}
{
}

void MyApp::drawEvent() {
    _plot.draw();
    swapBuffers();
};

void MyApp::mousePressEvent(MouseEvent& event) {
    if(event.button() != MouseEvent::Button::Left) return;
    _plot.setPosition(event.position());
    event.setAccepted();
}

void MyApp::mouseReleaseEvent(MouseEvent& event) {
    event.setAccepted();
}

void MyApp::mouseMoveEvent(MouseMoveEvent& event) {
    if(!(event.buttons() & MouseMoveEvent::Button::Left)) return;
    _plot.setMovePosition(event.position());
    event.setAccepted();
    redraw();
}

MAGNUM_APPLICATION_MAIN(MyApp)
