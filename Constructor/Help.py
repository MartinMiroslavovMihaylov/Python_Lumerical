import json
from pathlib import Path
import textwrap


class HelpDatabase:
    BASE = Path(__file__).parent / "help"

    # --- Initialize empty dicts first ---
    STRUCTURES = {}
    SOLVERS = {}
    STARSOLVER = {}
    RESULTS = {}
    FUNCTIONS = {}

    # --- Load JSON helper ---
    @staticmethod
    def load(name):
        path = HelpDatabase.BASE / name
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            with open(path, "w", encoding="utf-8") as f:
                json.dump({}, f)

        with open(path, encoding="utf-8") as f:
            return json.load(f)

    # --- Initialize after class exists ---
    @classmethod
    def initialize(cls):
        cls.STARTLUMERICAL = cls.load("start_lumerical.json")
        cls.STRUCTURES = cls.load("structures.json")
        cls.SOLVERS = cls.load("solvers.json")
        cls.STARSOLVER = cls.load("start_sim.json")
        cls.RESULTS = cls.load("results.json")
        cls.FUNCTIONS = cls.load("functions.json")

    # --- Add helpers ---
    @staticmethod
    def add_structure(idx, text):
        idx = str(idx)
        HelpDatabase.STRUCTURES[idx] = text
        path = HelpDatabase.BASE / "structures.json"
        with open(path, "w", encoding="utf-8") as f:
            json.dump(HelpDatabase.STRUCTURES, f, indent=2, ensure_ascii=False)

    @staticmethod
    def add_solver(idx, text):
        idx = str(idx)
        HelpDatabase.SOLVERS[idx] = text
        path = HelpDatabase.BASE / "solvers.json"
        with open(path, "w", encoding="utf-8") as f:
            json.dump(HelpDatabase.SOLVERS, f, indent=2, ensure_ascii=False)

    @staticmethod
    def add_results(idx, text):
        idx = str(idx)
        HelpDatabase.RESULTS[idx] = text
        path = HelpDatabase.BASE / "results.json"
        with open(path, "w", encoding="utf-8") as f:
            json.dump(HelpDatabase.RESULTS, f, indent=2, ensure_ascii=False)

    @staticmethod
    def add_functions(idx, text):
        idx = str(idx)
        HelpDatabase.FUNCTIONS[idx] = text
        path = HelpDatabase.BASE / "functions.json"
        with open(path, "w", encoding="utf-8") as f:
            json.dump(HelpDatabase.FUNCTIONS, f, indent=2, ensure_ascii=False)


# --- IMPORTANT: Load JSON files once class exists ---
HelpDatabase.initialize()


class HelpSystem:
    HELP_MAP = {
        "Start Lumerical": HelpDatabase.STARTLUMERICAL,
        "Structures": HelpDatabase.STRUCTURES,
        "Solvers": HelpDatabase.SOLVERS,
        "Start Simulation": HelpDatabase.STARSOLVER,
        "Results": HelpDatabase.RESULTS,
        "Functions": HelpDatabase.FUNCTIONS,
    }

    MENUS = {
        "Start Lumerical": {
            1: "Start Lumerical",
        },
        "Structures": {
            1: "Waveguide",
            2: "MMI2x1",
            3: "MMI2x2",
            4: "InverseTaper",
            5: "Directional Coupler",
            6: "StraightWaveguide",
            7: "WDM (Wavelength Division Multiplexing, angled MMI)",
            8: "BendWaveguide",
            9: "ArcWaveguide",
            10: "GratingCoupler",
            11: "RingGratingCoupler",
            12: "MZM",
            13: "Waveguide CHARGE",
        },
        "Solvers": {
            1: "Waveguide FDE",
            2: "MMI2x1 EME",
            3: "MMI2x2 EME",
            4: "InverseTaper EME",
            5: "Directional Coupler EME",
            6: "StraightWaveguide EME",
            7: "WDM angled MMI - EME",
            8: "MMI2x1 FDTD",
            9: "MMI2x2 FDTD",
            10: "InverseTaper FDTD",
            11: "Directional Coupler FDTD",
            12: "StraightWaveguide FDTD",
            13: "BendWaveguide FDTD",
            14: "ArcWaveguide FDTD",
            15: "GratingCoupler FDTD",
            16: "RingGratingCoupler FDTD",
            17: "MZM CHARGE",
            18: "MZM FEEM",
        },
        "Start Simulation": {
            1: "Start FDE Simulation",
            2: "Start EME Simulation",
            3: "Start FDTD Simulation",
            4: "Start CHARGE Simulation",
            5: "Start FEEM Simulation",
        },
        "Results": {
            1: "Extract FDE Results",
            2: "Extract and set Overlap Analysis with FDE Solver",
            3: "Extract EME Results",
            4: "Extract FDTD Results",
        },
        "Functions": {
            1: "Remove all objects and solvers from Lumerical workspace",
            2: "Terminal loading bar progress animation",
            3: "Create simulation log file from Lumerical setup",
        },
    }

    STRUCTURE_SOLVERS = {
        1: "FDE only",
        2: "EME & FDTD",
        3: "EME & FDTD",
        4: "EME & FDTD",
        5: "EME & FDTD",
        6: "EME & FDTD",
        7: "EME only",
        8: "FDTD only",
        9: "FDTD only",
        10: "FDTD only",
        11: "FDTD only",
        12: "CHARGE only",
        13: "CHARGE only",
    }

    @staticmethod
    def menu_help():
        print("=== Help Menu ===\n")

        major_categories = [
            "Start Lumerical",
            "Structures",
            "Solvers",
            "Start Simulation",
            "Results",
            "Functions",
        ]

        for idx, cat in enumerate(major_categories, 1):
            print(f"{idx}: {cat}")

        print("\nUsage hint:")
        print('HelpSystem.show_category("<CategoryName>", <ItemNumber>)')
        print('Example: HelpSystem.show_category("Structures", 1)')

    @staticmethod
    def show_category(category_name, item_number=None):
        if category_name not in HelpSystem.MENUS:
            print(f"Category '{category_name}' not found.")
            return

        menu_items = HelpSystem.MENUS[category_name]

        if item_number is None:
            print(f"\n--- {category_name} ---")
            for idx, name in menu_items.items():
                solver_info = ""
                if category_name == "Structures":
                    solver_info = HelpSystem.STRUCTURE_SOLVERS.get(idx, "")
                    if solver_info:
                        solver_info = f" - {solver_info}"
                print(f"{idx}: {name}{solver_info}")
        else:
            if item_number not in menu_items:
                print(f"Item {item_number} not found in {category_name}")
                return

            help_text = HelpSystem.HELP_MAP.get(category_name, {}).get(
                str(item_number),
                "No detailed help available.",
            )

            print(f"\n--- Detailed Help for {menu_items[item_number]} ---\n")
            print(textwrap.dedent(help_text))