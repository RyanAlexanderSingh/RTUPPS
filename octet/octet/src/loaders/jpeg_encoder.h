////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2012-2014
//
// Modular Framework for OpenGLES2 rendering on multiple platforms.
//
//
// jpeg file encoder - tiny and fast
//
// See http://en.wikipedia.org/wiki/JPEG
// 
namespace octet { namespace loaders {
  class jpeg_encoder {
  public:
    jpeg_encoder() {
    }

    bool encode(dynarray<uint8_t> data, uint32_t width, uint32_t height, int stride, const uint8_t *src) {
      /*if ((height & 7) != 0 || (width & 7) != 0) {
        return false;
      }

      // guess
      data.reserve(width * height / 2);

      static const uint8_t hdr[] = {
        0xff, 0xd8, // SOI
        0xff, 0xe0, 0x00, 0x10, // APP0
          'J', 'F', 'I', 'F', 0x00, 0x01, 0x01, 0x01, 0x00, 0x48, 0x00, 0x48, 0x00, 0x00,

        0xff, 0xdb, 0x00, 0x43, // DQT
          0x00,
          0x03, 0x02, 0x02, 0x03, 0x02, 0x02, 0x03, 0x03,
          0x03, 0x03, 0x04, 0x03, 0x03, 0x04, 0x05, 0x08,
          0x05, 0x05, 0x04, 0x04, 0x05, 0x0a, 0x07, 0x07,
          0x06, 0x08, 0x0c, 0x0a, 0x0c, 0x0c, 0x0b, 0x0a,
          0x0b, 0x0b, 0x0d, 0x0e, 0x12, 0x10, 0x0d, 0x0e,
          0x11, 0x0e, 0x0b, 0x0b, 0x10, 0x16, 0x10, 0x11,
          0x13, 0x14, 0x15, 0x15, 0x15, 0x0c, 0x0f, 0x17,
          0x18, 0x16, 0x14, 0x18, 0x12, 0x14, 0x15, 0x14,
        0xff, 0xdb, 0x00, 0x43, // DQT
          0x01,
          0x03, 0x04, 0x04, 0x05, 0x04, 0x05, 0x09, 0x05,
          0x05, 0x09, 0x14, 0x0d, 0x0b, 0x0d, 0x14, 0x14,
          0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
          0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
          0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
          0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
          0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
          0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14, 0x14,
        0xff, 0xc0, 0x00, 0x11, // SOF0
          0x08, 0x02, 0x00, 0x02, 0x00, 0x03, 0x01, 0x11, 0x00, 0x02, 0x11, 0x01, 0x03, 0x11, 0x01,
        0xff, 0xc4, 0x00, 0x1e, // DHT
          0x00, 0x01, 0x00, 0x03, 0x00, 0x02, 0x03, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
          0x00, 0x00, 0x01, 0x02, 0x03, 0x04, 0x06, 0x05, 0x07, 0x09, 0x08, 0x0a,
        0xff, 0xc4, 0x00, 0x55, // DHT
          0x10, 0x00, 0x02, 0x00, 0x03, 0x04, 0x05, 0x08, 0x06, 0x06, 0x04, 0x0b, 0x05, 0x07, 0x05, 0x00,
          0x00, 0x00, 0x01, 0x02, 0x03, 0x11, 0x04, 0x05, 0x21, 0x31, 0x06, 0x12, 0x41, 0x51, 0x61, 0x07,
          0x14, 0x53, 0x71, 0x81, 0x92, 0xb1, 0xf0, 0x08, 0x13, 0x34, 0x72, 0x91, 0xa1, 0x09, 0x22, 0x32,
          0x42, 0xc1, 0xd1, 0x23, 0x52, 0xe1, 0xf1, 0x15, 0x16, 0x17, 0x36, 0x43, 0x62, 0x73, 0x82, 0x94,
          0xb2, 0xd2, 0x24, 0x55, 0x63, 0xa3, 0xc2, 0x18, 0x33, 0x35, 0x64, 0x74, 0x93, 0xa2, 0x25, 0x44,
          0x45, 0x54, 0x83,
        0xff, 0xc4, 0x00, 0x1d, // DHT
          0x01, 0x01, 0x01, 0x01, 0x00, 0x02, 0x03, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
          0x00, 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09,
        0xff, 0xc4, 0x00, 0x43, // DHT
          0x11, 0x01, 0x00, 0x01, 0x02, 0x03, 0x05, 0x02, 0x0b, 0x06, 0x04, 0x05, 0x04, 0x03, 0x01, 0x00,
          0x00, 0x00, 0x01, 0x02, 0x11, 0x03, 0x04, 0x31, 0x05, 0x12, 0x21, 0x32, 0x51, 0x06, 0x41, 0x07,
          0x13, 0x14, 0x22, 0x61, 0x71, 0x81, 0x91, 0xa1, 0xb1, 0xd1, 0x15, 0x34, 0x42, 0x53, 0x72, 0xf0,
          0x52, 0x62, 0xc1, 0xe1, 0x23, 0x82, 0x92, 0xa2, 0xc2, 0x24, 0x43, 0xb2, 0xd2, 0x33, 0x63, 0x73,
          0x93,
        0xff, 0xda, 0x00, 0x0c, // SOS
          0x03, 0x01, 0x00, 0x02, 0x11, 0x03, 0x11, 0x00, 0x3f, 0x00, 0xfd, 0x7f, 
      };*/
      return 0;
    }
  };
}}

