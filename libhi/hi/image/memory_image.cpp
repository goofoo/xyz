#include <hi/image/memory_image.hpp>

namespace hi {
memory_image::memory_image()
    : width_(0),
      height_(0),
      pixels_(nullptr),
      hDc_(nullptr),
      hDib_(nullptr),
      hOldDib_(nullptr) {}

bool memory_image::create(hi::sint width, hi::sint height) {
  memory_image image;  // 一時オブジェクト
  if (image.init(width, height)) {
    this->swap(image);
    return true;
  } else {
    return false;
  }
}

void memory_image::destroy() {
  if (hOldDib_) {
    ::SelectObject(hDc_, hOldDib_);
    hOldDib_ = nullptr;
  }

  if (hDc_) {
    ::DeleteDC(hDc_);
    hDc_ = nullptr;
  }

  if (hDib_) {
    ::DeleteObject(hDib_);
    hDib_ = nullptr;
  }
}

void memory_image::clear() {
  std::memset(pixels_, 0, sizeof(hi::dword) * width_ * height_);
}

bool memory_image::load(std::tchar_t const *fileName) {
  memory_image temp;          // 一時ビットマップ
  bool is_completed = false;  // 正しくロードできたかどうかのフラグ
  HDC hdc = nullptr;
  HGDIOBJ old = nullptr;

  HANDLE hBitmap =
      ::LoadImage(nullptr, fileName, IMAGE_BITMAP, 0, 0, LR_LOADFROMFILE);
  if (!hBitmap) {
    goto __memory_image_load_exit;
  }

  BITMAP bitmap;
  ::GetObject(hBitmap, sizeof(bitmap), &bitmap);

  if (!temp.init(bitmap.bmWidth, bitmap.bmHeight)) {
    goto __memory_image_load_exit;
  }

  hdc = ::CreateCompatibleDC(temp.hDc_);
  if (!hdc) {
    goto __memory_image_load_exit;
  }

  old = ::SelectObject(hdc, hBitmap);
  if (!old) {
    goto __memory_image_load_exit;
  }

  is_completed = (0 != ::BitBlt(temp.hDc_, 0, 0, temp.width(), temp.height(),
                                hdc, 0, 0, SRCCOPY));

  if (is_completed) {
    this->swap(temp);
  }

__memory_image_load_exit:
  if (hdc) {
    if (old) {
      ::SelectObject(hdc, old);
    }
    ::DeleteDC(hdc);
  }
  if (hBitmap) {
    ::DeleteObject(hBitmap);
  }

  return is_completed;
}

bool memory_image::copy(memory_image const &image) {
  memory_image temp;
  if (!temp.init(image.width(), image.height())) {
    return false;
  }

  std::memcpy(temp.pixels_, image.pixels_,
              sizeof(hi::dword) * image.width() * image.height());
  this->swap(temp);

  return true;
}

void memory_image::raster(HDC hdc) const {
  ::BitBlt(hdc, 0, 0, width_, height_, hDc_, 0, 0, SRCCOPY);
}

void memory_image::swap(memory_image &image) {
  std::swap(width_, image.width_);
  std::swap(height_, image.height_);
  std::swap(hDib_, image.hDib_);
  std::swap(hDc_, image.hDc_);
  std::swap(hOldDib_, image.hOldDib_);
  std::swap(pixels_, image.pixels_);
}

bool memory_image::init(hi::sint width, hi::sint height) {
  width_ = width;
  height_ = height;

  BITMAPINFO bmi;
  std::memset(&bmi, 0, sizeof(bmi));
  bmi.bmiHeader.biSize = sizeof(bmi);
  bmi.bmiHeader.biWidth = width;
  bmi.bmiHeader.biHeight = -height;
  bmi.bmiHeader.biPlanes = 1;
  bmi.bmiHeader.biBitCount = 32;
  bmi.bmiHeader.biCompression = BI_RGB;

  hDib_ = ::CreateDIBSection(nullptr, &bmi, DIB_RGB_COLORS,
                             reinterpret_cast<void **>(&pixels_), nullptr, 0);
  if (!hDib_) {
    destroy();
    return false;
  }

  hDc_ = ::CreateCompatibleDC(nullptr);
  if (!hDc_) {
    destroy();
    return false;
  }

  hOldDib_ = ::SelectObject(hDc_, hDib_);
  if (!hOldDib_) {
    destroy();
    return false;
  }

  return true;
}

}  // end of namespace hi
