# Building libhex9 (and GeoPlegma) on Windows

If `cargo build` for GeoPlegma failed on Windows somewhere inside `hex9-sys`
with errors from `h9_warp_embedded.cpp` — something about `__asm__`, `.incbin`,
or a bare `C2059: syntax error` — you have hit a known limitation, and it is
not your fault. This page explains what is going on and gives you two ways
through, one of which takes about ten minutes.

## What's actually wrong

libhex9 embeds a ~binary data blob (its projection warp field) directly into
the library at compile time using the GNU assembler's `.incbin` directive.
That keeps the library self-contained — no data file to find at runtime — but
GNU inline assembly is something **MSVC does not support**.

Rust on Windows defaults to the MSVC toolchain (`x86_64-pc-windows-msvc`), and
`hex9-sys` builds libhex9's C++ core with whatever C++ compiler matches your
Rust target. Default Rust + default Visual Studio compiler → the `.incbin`
line cannot compile. Nothing about your setup is broken; the MSVC route simply
doesn't exist for this library.

Two working routes:

## Route A — WSL2 (recommended: fastest and most reliable)

Everything builds cleanly on Linux, and WSL2 gives you a real Ubuntu inside
Windows with almost no friction.

1. In an **administrator** PowerShell:

   ```powershell
   wsl --install -d Ubuntu
   ```

   Reboot if asked, then open the "Ubuntu" app.

2. Inside Ubuntu, install the toolchain (`libclang-dev` is for bindgen, which
   generates the Rust↔C bindings):

   ```sh
   sudo apt update
   sudo apt install -y build-essential cmake clang libclang-dev git
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   source "$HOME/.cargo/env"
   ```

3. Clone and build **inside the Linux filesystem** (your home directory, *not*
   `/mnt/c/...` — builds on `/mnt/c` are painfully slow):

   ```sh
   cd ~
   git clone <your GeoPlegma fork/repo>
   cd GeoPlegma
   cargo build
   ```

That's it. VS Code works beautifully with WSL2 if you install the "WSL"
extension — you edit in Windows, build in Linux, and never notice the seam.

## Route B — native Windows with the GNU toolchain

If you specifically need native Windows binaries, switch Rust to the GNU
target and use MinGW-w64 as the C++ compiler. This route is sound in
principle but much less travelled than Route A — if it fights you, don't
burn a day on it; take Route A.

1. Switch Rust to the GNU toolchain:

   ```powershell
   rustup toolchain install stable-x86_64-pc-windows-gnu
   rustup default stable-x86_64-pc-windows-gnu
   ```

   (Or, to switch only for this project, run
   `rustup override set stable-x86_64-pc-windows-gnu` inside the repo.)

2. Install [MSYS2](https://www.msys2.org/), then in the **UCRT64** shell:

   ```sh
   pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-cmake \
             mingw-w64-ucrt-x86_64-ninja
   ```

   Add `C:\msys64\ucrt64\bin` to your Windows `PATH` so cargo can find
   `g++`, `cmake`, and `ninja`.

3. bindgen needs libclang. Install LLVM and point at it:

   ```powershell
   winget install LLVM.LLVM
   setx LIBCLANG_PATH "C:\Program Files\LLVM\bin"
   ```

4. Tell CMake to use Ninja (avoids generator-guessing), then build:

   ```powershell
   setx CMAKE_GENERATOR Ninja
   # open a fresh terminal so the setx values are picked up
   cargo build
   ```

Keep the whole dependency tree on the `-gnu` target — mixing MSVC-built and
MinGW-built C++ libraries in one binary does not work.

## Sanity check: build libhex9 by itself

If GeoPlegma still won't build, check whether libhex9 alone does — it isolates
whether the problem is this library or something downstream:

```sh
git clone https://github.com/MrBenGriffin/libhex9
cd libhex9
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

All tests passing means libhex9 is fine on your machine, and the issue lives
in the GeoPlegma layer — in which case do get in touch with the error output,
which is genuinely welcome.
