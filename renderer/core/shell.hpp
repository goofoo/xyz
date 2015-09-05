#ifndef __TGIR_CORE_SHELL_HPP__
#define __TGIR_CORE_SHELL_HPP__

namespace tgir {
class Scene;

std::tstring const &GetOutdirPath();

int GetRenderingTime();
int GetIntervalTime();

bool IsFinished();
void FinishRendering();

void Display();
void Resize(int, int);
void Keyboard(unsigned char, int, int);

int Shell(std::tistream &in, bool const gui = true);

}  // end of namespace tgir

#endif __TGIR_CORE_SHELL_HPP__
