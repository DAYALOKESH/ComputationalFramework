"""
DEM Profile Extraction API
Fetches elevation data from Google Earth Engine along a straight line between two coordinates.
Returns distance-height pairs at equally spaced intervals. 

Requirements:
    pip install earthengine-api

Setup:
    1. Create a Google Cloud Project:  https://console.cloud.google.com/
    2. Enable the Earth Engine API for your project
    3. Run: python -c "import ee; ee. Authenticate()"
    4. Use your project ID when initializing
"""

import ee
import math
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict, Any
from enum import Enum
import json
import os


class DEMDataset(Enum):
    """Available DEM datasets in Google Earth Engine."""
    SRTM_30M = "USGS/SRTMGL1_003"
    SRTM_90M = "CGIAR/SRTM90_V4"
    ALOS_30M = "JAXA/ALOS/AW3D30/V3_2"
    COPERNICUS_30M = "COPERNICUS/DEM/GLO30"
    NASADEM = "NASA/NASADEM_HGT/001"


@dataclass
class Coordinate:
    """Represents a geographic coordinate."""
    latitude: float
    longitude: float
    
    def __post_init__(self):
        if not -90 <= self. latitude <= 90:
            raise ValueError(f"Latitude must be between -90 and 90, got {self.latitude}")
        if not -180 <= self. longitude <= 180:
            raise ValueError(f"Longitude must be between -180 and 180, got {self.longitude}")
    
    def to_list(self) -> List[float]:
        """Returns [longitude, latitude] for GEE compatibility."""
        return [self.longitude, self.latitude]
    
    def to_dict(self) -> Dict[str, float]:
        return {"latitude": self. latitude, "longitude":  self.longitude}


@dataclass
class ElevationPoint:
    """Represents a single elevation point along the profile."""
    distance_meters: float
    elevation_meters: Optional[float]
    latitude: float
    longitude: float
    
    def to_dict(self) -> Dict[str, Any]: 
        return {
            "distance_m": round(self.distance_meters, 2),
            "elevation_m": round(self.elevation_meters, 2) if self.elevation_meters is not None else None,
            "latitude": round(self.latitude, 6),
            "longitude": round(self. longitude, 6)
        }


@dataclass
class DEMProfileResult:
    """Result container for DEM profile extraction."""
    start_point:  Coordinate
    end_point: Coordinate
    total_distance_meters:  float
    num_points: int
    interval_meters: float
    dataset_used: str
    elevation_points: List[ElevationPoint] = field(default_factory=list)
    statistics: Dict[str, float] = field(default_factory=dict)
    success: bool = True
    error_message: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "success": self.success,
            "error_message": self. error_message,
            "metadata": {
                "start_point": self.start_point.to_dict(),
                "end_point": self. end_point.to_dict(),
                "total_distance_m": round(self.total_distance_meters, 2),
                "num_points": self.num_points,
                "interval_m": round(self. interval_meters, 2),
                "dataset":  self.dataset_used
            },
            "statistics": self.statistics,
            "profile": [p.to_dict() for p in self. elevation_points]
        }
    
    def to_json(self, indent:  int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)
    
    def get_distance_elevation_pairs(self) -> List[Tuple[float, Optional[float]]]: 
        """Returns list of (distance, elevation) tuples."""
        return [(p.distance_meters, p.elevation_meters) for p in self.elevation_points]
    
    def to_text(self) -> str:
        """Returns distance height pairs separated by space, one per line."""
        lines = []
        for p in self.elevation_points:
            dist = round(p.distance_meters, 2)
            elev = round(p.elevation_meters, 2) if p.elevation_meters is not None else "NaN"
            lines.append(str(dist) + " " + str(elev))
        return "\n".join(lines)
    
    def save(self, filepath: str) -> None:
        """Save distance height pairs to file."""
        with open(filepath, 'w') as f:
            f.write(self.to_text())
        print(f"Saved to {filepath}")


class DEMProfileExtractor:
    """
    Extracts elevation profiles from Google Earth Engine DEM data.
    """
    
    _global_initialized = False
    _current_project = None
    
    def __init__(
        self, 
        dataset:  DEMDataset = DEMDataset. SRTM_30M, 
        project:  Optional[str] = None,
        auto_init: bool = True
    ):
        self.dataset = dataset
        self.project = project or os.environ.get('EE_PROJECT_ID') or os.environ.get('GOOGLE_CLOUD_PROJECT')
        
        if auto_init: 
            self._initialize_ee()
    
    def _initialize_ee(self) -> None:
        """Initialize Google Earth Engine with project."""
        if DEMProfileExtractor._global_initialized and DEMProfileExtractor._current_project == self.project:
            return
        
        if not self.project:
            print("\n" + "="*70)
            print("GOOGLE EARTH ENGINE SETUP REQUIRED")
            print("="*70)
            print("\nTo use this tool, you need a Google Cloud Project with Earth Engine API enabled.")
            print("\nSetup steps:")
            print("1. Go to:  https://console.cloud.google.com/")
            print("2. Create a new project (or select existing)")
            print("3. Enable the 'Earth Engine API' for your project")
            print("4. Register your project at: https://code.earthengine. google.com/register")
            print("\nYou can also set the EE_PROJECT_ID environment variable to avoid this prompt.")
            print("="*70 + "\n")
            
            self.project = input("Enter your Google Cloud Project ID: ").strip()
            
            if not self. project:
                raise ValueError("Project ID is required.")
        
        try:
            ee. Initialize(project=self.project)
            DEMProfileExtractor._global_initialized = True
            DEMProfileExtractor._current_project = self.project
            print(f"Earth Engine initialized with project: {self.project}")
            
        except ee.EEException as e:
            error_msg = str(e)
            
            if "authenticate" in error_msg.lower() or "credentials" in error_msg.lower():
                print("Authentication required. Running authentication flow...")
                try:
                    ee.Authenticate()
                    ee.Initialize(project=self.project)
                    DEMProfileExtractor._global_initialized = True
                    DEMProfileExtractor._current_project = self.project
                    print(f"Earth Engine initialized with project: {self. project}")
                except Exception as auth_error:
                    raise RuntimeError(f"Authentication failed: {auth_error}")
            else:
                raise RuntimeError(f"Failed to initialize Earth Engine:  {e}")
    
    @staticmethod
    def haversine_distance(coord1: Coordinate, coord2: Coordinate) -> float:
        """Calculate the great-circle distance between two points in meters."""
        R = 6371000
        
        lat1_rad = math.radians(coord1.latitude)
        lat2_rad = math.radians(coord2.latitude)
        delta_lat = math.radians(coord2.latitude - coord1.latitude)
        delta_lon = math.radians(coord2.longitude - coord1.longitude)
        
        a = (math.sin(delta_lat / 2) ** 2 + 
             math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(delta_lon / 2) ** 2)
        c = 2 * math. atan2(math. sqrt(a), math.sqrt(1 - a))
        
        return R * c
    
    def _interpolate_great_circle(
        self,
        start:  Coordinate,
        end: Coordinate,
        num_points: int
    ) -> List[Coordinate]:
        """Generate equally spaced coordinates along a great circle path."""
        if num_points < 2:
            raise ValueError("num_points must be at least 2")
        
        lat1 = math.radians(start.latitude)
        lon1 = math. radians(start. longitude)
        lat2 = math. radians(end. latitude)
        lon2 = math. radians(end. longitude)
        
        d = 2 * math.asin(math.sqrt(
            math.sin((lat2 - lat1) / 2) ** 2 +
            math.cos(lat1) * math.cos(lat2) * math.sin((lon2 - lon1) / 2) ** 2
        ))
        
        coordinates = []
        
        for i in range(num_points):
            fraction = i / (num_points - 1) if num_points > 1 else 0
            
            if d < 1e-10: 
                coordinates.append(Coordinate(latitude=start.latitude, longitude=start.longitude))
                continue
            
            A = math.sin((1 - fraction) * d) / math.sin(d)
            B = math.sin(fraction * d) / math.sin(d)
            
            x = A * math. cos(lat1) * math.cos(lon1) + B * math.cos(lat2) * math.cos(lon2)
            y = A * math.cos(lat1) * math.sin(lon1) + B * math.cos(lat2) * math.sin(lon2)
            z = A * math.sin(lat1) + B * math.sin(lat2)
            
            lat = math.degrees(math.atan2(z, math.sqrt(x ** 2 + y ** 2)))
            lon = math. degrees(math. atan2(y, x))
            
            coordinates.append(Coordinate(latitude=lat, longitude=lon))
        
        return coordinates
    
    def _get_elevation_band(self) -> str:
        """Get the elevation band name for the current dataset."""
        band_mapping = {
            DEMDataset. SRTM_30M: "elevation",
            DEMDataset.SRTM_90M:  "elevation",
            DEMDataset. ALOS_30M: "DSM",
            DEMDataset. COPERNICUS_30M: "DEM",
            DEMDataset.NASADEM: "elevation"
        }
        return band_mapping.get(self. dataset, "elevation")
    
    def _fetch_elevations_batch(self, coordinates:  List[Coordinate], batch_size: int = 100) -> List[Optional[float]]: 
        """Fetch elevation values for coordinates in batches."""
        dem = ee.Image(self.dataset.value)
        band = self._get_elevation_band()
        
        all_elevations:  List[Optional[float]] = [None] * len(coordinates)
        
        for batch_start in range(0, len(coordinates), batch_size):
            batch_end = min(batch_start + batch_size, len(coordinates))
            batch_coords = coordinates[batch_start: batch_end]
            
            points = ee. FeatureCollection([
                ee.Feature(ee.Geometry.Point(coord.to_list()), {'index': i})
                for i, coord in enumerate(batch_coords)
            ])
            
            sampled = dem.select(band).sampleRegions(
                collection=points,
                scale=30,
                geometries=False
            )
            
            try:
                results = sampled.getInfo()
                
                for feature in results. get('features', []):
                    props = feature.get('properties', {})
                    idx = props. get('index')
                    elev = props.get(band)
                    if idx is not None and elev is not None: 
                        all_elevations[batch_start + idx] = float(elev)
                        
            except Exception as e:
                print(f"Warning: Failed to fetch batch {batch_start}-{batch_end}:  {e}")
        
        return all_elevations
    
    def _calculate_statistics(self, elevations: List[Optional[float]]) -> Dict[str, Any]:
        """Calculate statistics for the elevation profile."""
        valid_elevations = [e for e in elevations if e is not None]
        
        if not valid_elevations:
            return {"error": "No valid elevation data"}
        
        min_elev = min(valid_elevations)
        max_elev = max(valid_elevations)
        mean_elev = sum(valid_elevations) / len(valid_elevations)
        
        total_ascent = 0.0
        total_descent = 0.0
        prev_elev = None
        
        for elev in elevations:
            if elev is not None and prev_elev is not None:
                diff = elev - prev_elev
                if diff > 0:
                    total_ascent += diff
                else:
                    total_descent += abs(diff)
            if elev is not None:
                prev_elev = elev
        
        return {
            "min_elevation_m": round(min_elev, 2),
            "max_elevation_m": round(max_elev, 2),
            "mean_elevation_m": round(mean_elev, 2),
            "elevation_range_m": round(max_elev - min_elev, 2),
            "total_ascent_m": round(total_ascent, 2),
            "total_descent_m": round(total_descent, 2),
            "valid_points": len(valid_elevations),
            "invalid_points": len(elevations) - len(valid_elevations)
        }
    
    def get_profile(
        self,
        start: Tuple[float, float],
        end:  Tuple[float, float],
        num_points: int = 100
    ) -> DEMProfileResult:
        """
        Extract elevation profile between two points.
        
        Args: 
            start: Starting point as (latitude, longitude)
            end:  Ending point as (latitude, longitude)
            num_points: Number of equally spaced points to sample
            
        Returns:
            DEMProfileResult containing the elevation profile
        """
        try:
            start_coord = Coordinate(latitude=start[0], longitude=start[1])
            end_coord = Coordinate(latitude=end[0], longitude=end[1])
            
            total_distance = self.haversine_distance(start_coord, end_coord)
            interval = total_distance / (num_points - 1) if num_points > 1 else 0
            
            print(f"Extracting profile: {total_distance/1000:.2f} km, {num_points} points...")
            
            coordinates = self._interpolate_great_circle(start_coord, end_coord, num_points)
            elevations = self._fetch_elevations_batch(coordinates)
            
            elevation_points = []
            for i, (coord, elev) in enumerate(zip(coordinates, elevations)):
                distance = interval * i
                elevation_points.append(ElevationPoint(
                    distance_meters=distance,
                    elevation_meters=elev,
                    latitude=coord.latitude,
                    longitude=coord. longitude
                ))
            
            statistics = self._calculate_statistics(elevations)
            
            print(f"Profile extracted.  Valid points: {statistics. get('valid_points', 0)}/{num_points}")
            
            return DEMProfileResult(
                start_point=start_coord,
                end_point=end_coord,
                total_distance_meters=total_distance,
                num_points=num_points,
                interval_meters=interval,
                dataset_used=self. dataset. value,
                elevation_points=elevation_points,
                statistics=statistics,
                success=True
            )
            
        except Exception as e: 
            return DEMProfileResult(
                start_point=Coordinate(start[0], start[1]),
                end_point=Coordinate(end[0], end[1]),
                total_distance_meters=0,
                num_points=num_points,
                interval_meters=0,
                dataset_used=self.dataset.value,
                success=False,
                error_message=str(e)
            )
    
    def get_profile_by_interval(
        self,
        start: Tuple[float, float],
        end:  Tuple[float, float],
        interval_meters: float = 100
    ) -> DEMProfileResult: 
        """Extract elevation profile with a specified sampling interval in meters."""
        start_coord = Coordinate(latitude=start[0], longitude=start[1])
        end_coord = Coordinate(latitude=end[0], longitude=end[1])
        
        total_distance = self.haversine_distance(start_coord, end_coord)
        num_points = max(2, int(total_distance / interval_meters) + 1)
        
        return self.get_profile(start, end, num_points)


def get_elevation_profile(
    start:  Tuple[float, float],
    end: Tuple[float, float],
    num_points: int = 100,
    project: Optional[str] = None,
    dataset: DEMDataset = DEMDataset. SRTM_30M
) -> DEMProfileResult: 
    """
    Quick function to get elevation profile between two points.
    
    Args: 
        start: Starting point as (latitude, longitude)
        end: Ending point as (latitude, longitude)
        num_points:  Number of sample points
        project:  Google Cloud Project ID
        dataset: DEM dataset to use
        
    Returns: 
        DEMProfileResult containing the elevation profile
    """
    extractor = DEMProfileExtractor(dataset=dataset, project=project)
    return extractor.get_profile(start, end, num_points)


if __name__ == "__main__":
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(description="Extract DEM elevation profile")
    parser.add_argument("--project", type=str, help="Google Cloud Project ID")
    parser.add_argument("--start-lat", type=float, required=True, help="Starting latitude")
    parser.add_argument("--start-lon", type=float, required=True, help="Starting longitude")
    parser.add_argument("--end-lat", type=float, required=True, help="Ending latitude")
    parser.add_argument("--end-lon", type=float, required=True, help="Ending longitude")
    parser.add_argument("--num-points", type=int, default=100, help="Number of sample points")
    parser.add_argument("--output", type=str, help="Output file path")
    
    args = parser.parse_args()
    
    extractor = DEMProfileExtractor(project=args.project)
    result = extractor.get_profile(
        start=(args.start_lat, args.start_lon),
        end=(args.end_lat, args.end_lon),
        num_points=args.num_points
    )
    
    if result.success:
        output = result.to_text()
        if args.output:
            result.save(args. output)
        else:
            print(output)
    else:
        print(f"Error: {result.error_message}")
        sys.exit(1)