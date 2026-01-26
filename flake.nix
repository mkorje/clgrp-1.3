{
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
  };

  outputs =
    {
      self,
      nixpkgs,
    }:
    let
      system = "x86_64-linux";
      pkgs = nixpkgs.legacyPackages.${system};
    in
    {
      packages.${system} = rec {
        optarith = pkgs.callPackage ./optarith.nix { primes = 10000; };
        qform = pkgs.callPackage ./qform.nix { inherit optarith; };
      };

      devShells.${system}.default = pkgs.mkShell {
        nativeBuildInputs = with pkgs; [
          autoconf
          automake
          libtool
          pkg-config
        ];

        buildInputs = with pkgs; [
          openmpi
          gmp
          pari
          self.packages.${system}.optarith
          self.packages.${system}.qform
        ];

        shellHook = ''
          cat << EOF
          make distclean
          autoreconf
          ./configure CC=mpicc CPPFLAGS="-DDEBUG -DKEEP_FILES -DWITH_PARI" LIBS="-lpari -lm"
          ./configure CC=mpicc CPPFLAGS="-DKEEP_FILES -DWITH_PARI" LIBS="-lpari -lm"
          ./configure CC=mpicc CPPFLAGS="-DKEEP_FILE" LIBS="-lm"
          make
          EOF
        '';
      };
    };
}
