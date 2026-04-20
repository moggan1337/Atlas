#include "vec3.hpp"

// Explicit template instantiations for common types
namespace atlas {
namespace math {

template Vec3 Vec3::lerp<>(
    const Vec3& a, const Vec3& b, float t
);

} // namespace math
} // namespace atlas
