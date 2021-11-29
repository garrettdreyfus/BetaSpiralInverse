with import <nixpkgs> {}; (pkgs.buildFHSUserEnv { name = "fhs";
	profile = ''
	      set -e
	      eval "$(micromamba shell hook -s bash)"
	      micromamba activate betaspiral
	      set +e

	    '';
    targetPkgs = pkgs: (builtins.concatLists [ [ micromamba ] [] [pkgs.which]]);
 }).env

