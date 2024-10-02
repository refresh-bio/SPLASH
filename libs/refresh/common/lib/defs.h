#if defined(_MSC_VER)  /* Visual Studio */
#define REFRESH_FORCE_INLINE __forceinline
#define REFRESH_NO_INLINE __declspec(noinline)
#elif defined(__GNUC__)
#define REFRESH_FORCE_INLINE __inline__ __attribute__((always_inline, unused))
#define REFRESH_NO_INLINE __attribute__((noinline))
#else
#define REFRESH_FORCE_INLINE
#define REFRESH_NO_INLINE
#endif
